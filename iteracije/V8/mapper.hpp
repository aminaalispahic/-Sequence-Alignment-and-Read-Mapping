// mapper.hpp: k-mer indeks, seeding, extension (Smith-Waterman) i izbor najboljeg mapiranja
#pragma once
#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include "align.hpp"
#include "fasta.hpp"

namespace sw {

inline int enc(char c) {
    const Base base = encodeBase(c);
    return isCanonical(base) ? base : -1;
}
inline std::string revComp(const std::string& s) {
    std::string r(s.rbegin(), s.rend());
    for (char& c : r) switch (c) { case 'A': c = 'T'; break; case 'C': c = 'G'; break;
                                   case 'G': c = 'C'; break; case 'T': c = 'A'; break; default: break; }
    return r;
}

// Oba lanca svakog reada koriste isti par granica u dva kontinuirana bafera.
struct ReadBatch {
    std::string bases, reverseBases;
    std::vector<size_t> offsets{0};
    std::vector<std::string> ids;

    size_t size() const { return ids.size(); }
    EncodedView sequence(size_t i, bool reverse = false) const {
        const std::string& buffer = reverse ? reverseBases : bases;
        return EncodedView{buffer}.substr(offsets[i], offsets[i + 1] - offsets[i]);
    }
};

inline ReadBatch makeBatch(std::vector<Record> records) {
    ReadBatch batch;
    size_t total = 0;
    for (const auto& record : records) total += record.seq.size();
    batch.bases.reserve(total);
    batch.reverseBases.reserve(total);
    batch.offsets.reserve(records.size() + 1);
    batch.ids.reserve(records.size());
    for (auto& record : records) {
        const size_t start = batch.bases.size();
        appendEncoded(record.seq, batch.bases);
        appendReverseComplement(EncodedView{batch.bases}.substr(start, record.seq.size()),
                                batch.reverseBases);
        batch.offsets.push_back(batch.bases.size());
        batch.ids.push_back(std::move(record.id));
    }
    return batch;
}

struct KmerPositions {
    const uint32_t* direct = nullptr;
    const uint64_t* packed = nullptr;
    std::ptrdiff_t count = 0;

    void prefetchFirst() const {
#if defined(__GNUC__) || defined(__clang__)
        if (direct) __builtin_prefetch(direct, 0, 1);
#endif
    }

    template <class Visit>
    void forEach(Visit visit) const {
        if (direct) {
            for (std::ptrdiff_t i = 0; i < count; ++i) visit(direct[i]);
        } else {
            for (std::ptrdiff_t i = 0; i < count; ++i) visit(static_cast<uint32_t>(packed[i]));
        }
    }
};

// Direktni indeks: positions[offsets[code] .. offsets[code + 1]) za svaki k-mer.
// Za k > 13 ostaje sortirani indeks zbog rasta tabele s 4^k.
struct KmerIndex {
    int k = 12;

    void build(EncodedView ref, int k_) {
        k = k_;
        if (k > 13 || ref.size() > std::numeric_limits<uint32_t>::max()) {
            offsets = std::vector<uint32_t>{};
            positions = std::vector<uint32_t>{};
            packed.clear();
            forEachKmer(ref, [&](uint64_t code, uint32_t pos) {
                packed.push_back((code << 32) | pos);
            });
            std::sort(packed.begin(), packed.end());
            return;
        }

        packed = std::vector<uint64_t>{};
        const size_t codeCount = size_t{1} << (2 * k);
        offsets.assign(codeCount + 1, 0);
        forEachKmer(ref, [&](uint64_t code, uint32_t) { ++offsets[code + 1]; });
        for (size_t i = 1; i < offsets.size(); ++i) offsets[i] += offsets[i - 1];

        positions.resize(offsets.back());
        forEachKmer(ref, [&](uint64_t code, uint32_t pos) {
            positions[offsets[code]++] = pos;
        });
        // Kursori sada pokazuju krajeve raspona. Pomjeri ih za jedno mjesto.
        std::move_backward(offsets.begin(), offsets.end() - 1, offsets.end());
        offsets.front() = 0;
    }

    void build(const std::string& ref, int k_) {
        const auto encoded = encodeSequence(ref);
        build(EncodedView{encoded}, k_);
    }

    size_t size() const { return positions.size() + packed.size(); }

    void prefetchOffsets(uint64_t code) const {
#if defined(__GNUC__) || defined(__clang__)
        if (!offsets.empty() && code < offsets.size() - 1)
            __builtin_prefetch(offsets.data() + code, 0, 1);
#else
        (void)code;
#endif
    }

    KmerPositions lookup(uint64_t code) const {
        if (offsets.empty()) {
            auto lo = std::lower_bound(packed.begin(), packed.end(), code << 32);
            auto hi = std::lower_bound(lo, packed.end(), (code + 1) << 32);
            if (lo == hi) return {};
            return {nullptr, &*lo, hi - lo};
        }
        if (code >= offsets.size() - 1) return {};
        const uint32_t begin = offsets[code];
        const uint32_t end = offsets[code + 1];
        if (begin == end) return {};
        return {positions.data() + begin, nullptr, static_cast<std::ptrdiff_t>(end - begin)};
    }

private:
    std::vector<uint32_t> offsets, positions;
    std::vector<uint64_t> packed;

    template <class Visit>
    void forEachKmer(EncodedView ref, Visit visit) const {
        const uint64_t mask = (1ULL << (2 * k)) - 1;
        uint64_t code = 0; int valid = 0;
        for (size_t i = 0; i < ref.size(); ++i) {
            const Base c = ref[i];
            if (!isCanonical(c)) { valid = 0; code = 0; continue; }
            code = ((code << 2) | (uint64_t)c) & mask;
            if (++valid >= k) visit(code, static_cast<uint32_t>(i - k + 1));
        }
    }
};

struct MapParams {
    int k = 12;
    int margin = 20;            // margina oko kandidata; ujedno tolerancija za pomak dijagonale (indeli)
    int minVotes = 2;           // minimalan broj seed pogodaka za kandidata
    int maxOcc = 500;           // k-mere s više pojavljivanja se preskaču (ponavljajuće regije)
    int maxCand = 5;            // najviše kandidata po lancu
    int stride = 1;             // razmak između seedova u readu
    double minScoreFrac = 0.5;  // minimalni skor = frac * dužina * match
    Scoring sc;
    int band = 12;               // polu-sirina oko seed dijagonale, -1 = puna matrica
    double relVotes = 0.5;       // udio najboljeg broja glasova na oba lanca
};

struct MapResult {
    std::string status = "unmapped";   // mapped | ambiguous | unmapped
    char strand = '*';
    int64_t pos = -1;                  // 0-indeksirana pozicija na referenci
    int score = 0;
    std::string cigar = "*";
    int candidates = 0;
};

using SeedCandidate = std::pair<int64_t, int>;

struct SeedWorkspace {
    struct Seed { uint64_t code; int64_t readOffset; };
    struct Match { KmerPositions positions; int64_t readOffset; };
    struct Cluster { int64_t diag; int votes; int bestRun; };

    std::vector<Seed> seeds;
    std::vector<Match> matches;
    std::vector<int64_t> diags;
    std::vector<Cluster> clusters;
    std::array<std::vector<SeedCandidate>, 2> candidates;
};

struct Stats {
    double seedSec = 0, extSec = 0;
    uint64_t cells = 0, extensions = 0, tracebacks = 0;
    SwWorkspace swWorkspace;
    uint64_t candidatesBeforeRelVotes = 0, candidatesAfterRelVotes = 0;
    SeedWorkspace seedWorkspace;
};

inline void findCandidates(EncodedView seq, const KmerIndex& idx, const MapParams& p,
                           SeedWorkspace& workspace, std::vector<SeedCandidate>& out) {
    auto& seeds = workspace.seeds;
    auto& matches = workspace.matches;
    auto& diags = workspace.diags;
    auto& cl = workspace.clusters;
    seeds.clear();
    matches.clear();
    diags.clear();
    cl.clear();
    out.clear();

    // Prvo izracunaj kodove k-mera i zatrazi njihove offsete.
    const int k = idx.k;
    const uint64_t mask = (1ULL << (2 * k)) - 1;
    uint64_t code = 0; int valid = 0;
    for (size_t i = 0; i < seq.size(); ++i) {
        const Base c = seq[i];
        if (!isCanonical(c)) { valid = 0; code = 0; continue; }
        code = ((code << 2) | (uint64_t)c) & mask;
        if (++valid < k) continue;
        const int64_t off = (int64_t)i - k + 1;
        if (off % p.stride) continue;
        seeds.push_back({code, off});
        idx.prefetchOffsets(code);
    }

    // Zatim procitaj raspone i zatrazi pocetak niza pozicija.
    for (const auto& seed : seeds) {
        const auto positions = idx.lookup(seed.code);
        if (positions.count == 0 || positions.count > p.maxOcc) continue;
        positions.prefetchFirst();
        matches.push_back({positions, seed.readOffset});
    }

    // Na kraju procitaj pozicije i izracunaj dijagonale.
    for (const auto& match : matches) {
        match.positions.forEach([&](uint32_t pos) {
            diags.push_back((int64_t)pos - match.readOffset);
        });
    }
    std::sort(diags.begin(), diags.end());
    int64_t last = 0;
    for (size_t i = 0; i < diags.size();) {
        size_t j = i; while (j < diags.size() && diags[j] == diags[i]) ++j;
        int cnt = (int)(j - i); int64_t d = diags[i];
        if (cl.empty() || d - last > p.margin) cl.push_back({d, cnt, cnt});
        else { auto& c = cl.back(); c.votes += cnt; if (cnt > c.bestRun) { c.bestRun = cnt; c.diag = d; } }
        last = d; i = j;
    }
    for (auto& c : cl) if (c.votes >= p.minVotes) out.push_back({c.diag, c.votes});
    std::sort(out.begin(), out.end(), [](auto& x, auto& y) { return x.second > y.second; });
    if ((int)out.size() > p.maxCand) out.resize(p.maxCand);
}

inline std::vector<SeedCandidate> findCandidates(EncodedView seq, const KmerIndex& idx, const MapParams& p) {
    SeedWorkspace workspace;
    std::vector<SeedCandidate> out;
    findCandidates(seq, idx, p, workspace, out);
    return out;
}

inline std::vector<std::pair<int64_t, int>> findCandidates(const std::string& seq, const KmerIndex& idx, const MapParams& p) {
    const auto encoded = encodeSequence(seq);
    return findCandidates(EncodedView{encoded}, idx, p);
}

struct CandidateWindow {
    EncodedView seq, region;
    int64_t refStart, diagonal;
    char strand;
    int score = 0;

    uint64_t cells() const { return (uint64_t)seq.size() * region.size(); }
};

inline int scoreCandidate(const CandidateWindow& candidate, const MapParams& p, Stats& st) {
    ++st.extensions;
    if (p.band >= 0)
        return smithWatermanScoreBanded(candidate.seq, candidate.region, candidate.diagonal,
                                       p.band, p.sc, st.swWorkspace, st.cells);
    st.cells += candidate.cells();
    return smithWatermanScore(candidate.seq, candidate.region, p.sc, st.swWorkspace);
}

inline Alignment tracebackCandidate(const CandidateWindow& candidate, const MapParams& p, Stats& st) {
    ++st.extensions;
    ++st.tracebacks;
    if (p.band >= 0)
        return smithWatermanBanded(candidate.seq, candidate.region, candidate.diagonal,
                                  p.band, p.sc, st.swWorkspace, st.cells);
    st.cells += candidate.cells();
    return smithWaterman(candidate.seq, candidate.region, p.sc, st.swWorkspace);
}

inline MapResult selectBestCandidate(std::vector<CandidateWindow>& windows, int ncand,
                                     size_t readLength, const MapParams& p, Stats& st) {
    MapResult result;
    result.candidates = ncand;
    if (windows.empty()) return result;

    const int minScore = (int)(p.minScoreFrac * readLength * p.sc.match);
    if (windows.size() > 1) {
        for (auto& window : windows) window.score = scoreCandidate(window, p, st);
        windows.erase(std::remove_if(windows.begin(), windows.end(),
            [](const CandidateWindow& window) { return window.score <= 0; }), windows.end());
        if (windows.empty()) return result;

        // Sortiraj sve pozitivne skorove istim redoslijedom i kriterijem kao V4.
        std::sort(windows.begin(), windows.end(),
            [](const CandidateWindow& x, const CandidateWindow& y) { return x.score > y.score; });
        if (windows.front().score < minScore) return result;
    }

    const auto& best = windows.front();
    const Alignment alignment = tracebackCandidate(best, p, st);
    if (alignment.score <= 0 || alignment.score < minScore) return result;

    result.status = "mapped";
    result.strand = best.strand;
    result.pos = best.refStart + alignment.startY;
    result.score = alignment.score;
    result.cigar.clear();
    if (alignment.startX > 0) result.cigar += std::to_string(alignment.startX) + "S";
    result.cigar += alignment.cigar;
    if ((int)best.seq.size() > alignment.endX)
        result.cigar += std::to_string((int)best.seq.size() - alignment.endX) + "S";

    for (size_t i = 1; i < windows.size() && windows[i].score == result.score; ++i) {
        const auto& tied = windows[i];
        if (tied.strand != result.strand) {
            result.status = "ambiguous";
            break;
        }
        const Alignment alternative = tracebackCandidate(tied, p, st);
        const int64_t position = tied.refStart + alternative.startY;
        if (std::llabs(position - result.pos) > p.margin) {
            result.status = "ambiguous";
            break;
        }
    }
    return result;
}

inline MapResult mapRead(EncodedView read, EncodedView rc, const KmerIndex& idx, EncodedView ref,
                         const MapParams& p, Stats& st) {
    using clk = std::chrono::steady_clock;
    auto secs = [](clk::time_point a, clk::time_point b) { return std::chrono::duration<double>(b - a).count(); };

    const auto seedStart = clk::now();
    auto& workspace = st.seedWorkspace;
    auto& candidates = workspace.candidates;
    findCandidates(read, idx, p, workspace, candidates[0]);
    findCandidates(rc, idx, p, workspace, candidates[1]);
    // Svaki niz vec je sortiran po glasovima i ogranicen na maxCand.
    const int forwardVotes = candidates[0].empty() ? 0 : candidates[0].front().second;
    const int reverseVotes = candidates[1].empty() ? 0 : candidates[1].front().second;
    const int bestVotes = std::max(forwardVotes, reverseVotes);
    const int threshold = std::max(p.minVotes, static_cast<int>(p.relVotes * bestVotes));
    st.candidatesBeforeRelVotes += candidates[0].size() + candidates[1].size();
    const auto extensionStart = clk::now();
    st.seedSec += secs(seedStart, extensionStart);

    std::vector<CandidateWindow> windows;
    int ncand = 0;
    for (int s = 0; s < 2; ++s) {                       // oba lanca: read i njegov reverse-complement
        const auto seq = (s == 0) ? read : rc;
        for (const auto& c : candidates[s]) {
            if (c.second < threshold) break;
            ++ncand;
            int64_t rs = std::max<int64_t>(0, c.first - p.margin);
            int64_t re = std::min<int64_t>((int64_t)ref.size(), c.first + (int64_t)seq.size() + p.margin);
            if (re <= rs) continue;
            const auto region = ref.substr((size_t)rs, (size_t)(re - rs));
            windows.push_back({seq, region, rs, c.first - rs, s == 0 ? '+' : '-'});
        }
    }
    st.candidatesAfterRelVotes += ncand;
    auto result = selectBestCandidate(windows, ncand, read.size(), p, st);
    st.extSec += secs(extensionStart, clk::now());
    return result;
}

inline MapResult mapRead(const std::string& read, const KmerIndex& idx, const std::string& ref,
                         const MapParams& p, Stats& st) {
    const auto encodedRead = encodeSequence(read);
    const auto encodedRef = encodeSequence(ref);
    std::string reverse;
    reverse.reserve(encodedRead.size());
    appendReverseComplement(EncodedView{encodedRead}, reverse);
    return mapRead(EncodedView{encodedRead}, EncodedView{reverse}, idx, EncodedView{encodedRef}, p, st);
}

}  // namespace sw
