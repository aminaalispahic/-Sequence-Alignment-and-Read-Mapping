// mapper.hpp: k-mer indeks, seeding, extension (Smith-Waterman) i izbor najboljeg mapiranja
#pragma once
#include <chrono>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>
#include "align.hpp"

namespace sw {

inline int enc(char c) {
    switch (c) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; default: return -1; }
}
inline std::string revComp(const std::string& s) {
    std::string r(s.rbegin(), s.rend());
    for (char& c : r) switch (c) { case 'A': c = 'T'; break; case 'C': c = 'G'; break;
                                   case 'G': c = 'C'; break; case 'T': c = 'A'; break; default: break; }
    return r;
}

// Indeks: sortiran niz (kmer << 32 | pozicija). Podržava k <= 16 i reference do 4 Gb.
struct KmerIndex {
    int k = 12;
    std::vector<uint64_t> e;
    void build(const std::string& ref, int k_) {
        k = k_; e.clear();
        const uint64_t mask = (1ULL << (2 * k)) - 1;
        uint64_t code = 0; int valid = 0;
        for (size_t i = 0; i < ref.size(); ++i) {
            int c = enc(ref[i]);
            if (c < 0) { valid = 0; code = 0; continue; }       // k-mere s N se preskaču
            code = ((code << 2) | (uint64_t)c) & mask;
            if (++valid >= k) e.push_back((code << 32) | (uint32_t)(i - k + 1));
        }
        std::sort(e.begin(), e.end());
    }
    std::pair<std::vector<uint64_t>::const_iterator, std::vector<uint64_t>::const_iterator>
    lookup(uint64_t code) const {
        auto lo = std::lower_bound(e.begin(), e.end(), code << 32);
        auto hi = std::lower_bound(lo, e.end(), (code + 1) << 32);
        return {lo, hi};
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
};

struct MapResult {
    std::string status = "unmapped";   // mapped | ambiguous | unmapped
    char strand = '*';
    int64_t pos = -1;                  // 0-indeksirana pozicija na referenci
    int score = 0;
    std::string cigar = "*";
    int candidates = 0;
};

struct Stats { double seedSec = 0, extSec = 0; uint64_t cells = 0, extensions = 0; };

// Seeding: sve k-mere reada -> pozicije u indeksu -> glasanje po dijagonali (p - offset).
inline std::vector<std::pair<int64_t, int>> findCandidates(const std::string& seq, const KmerIndex& idx, const MapParams& p) {
    std::vector<int64_t> diags;
    const int k = idx.k;
    const uint64_t mask = (1ULL << (2 * k)) - 1;
    uint64_t code = 0; int valid = 0;
    for (size_t i = 0; i < seq.size(); ++i) {
        int c = enc(seq[i]);
        if (c < 0) { valid = 0; code = 0; continue; }
        code = ((code << 2) | (uint64_t)c) & mask;
        if (++valid < k) continue;
        int64_t off = (int64_t)i - k + 1;
        if (off % p.stride) continue;
        auto r = idx.lookup(code);
        if (r.second - r.first > p.maxOcc) continue;
        for (auto it = r.first; it != r.second; ++it) diags.push_back((int64_t)(uint32_t)(*it) - off);
    }
    std::sort(diags.begin(), diags.end());
    struct Cl { int64_t diag; int votes; int bestRun; };
    std::vector<Cl> cl; int64_t last = 0;
    for (size_t i = 0; i < diags.size();) {
        size_t j = i; while (j < diags.size() && diags[j] == diags[i]) ++j;
        int cnt = (int)(j - i); int64_t d = diags[i];
        if (cl.empty() || d - last > p.margin) cl.push_back({d, cnt, cnt});
        else { Cl& c = cl.back(); c.votes += cnt; if (cnt > c.bestRun) { c.bestRun = cnt; c.diag = d; } }
        last = d; i = j;
    }
    std::vector<std::pair<int64_t, int>> out;
    for (auto& c : cl) if (c.votes >= p.minVotes) out.push_back({c.diag, c.votes});
    std::sort(out.begin(), out.end(), [](auto& x, auto& y) { return x.second > y.second; });
    if ((int)out.size() > p.maxCand) out.resize(p.maxCand);
    return out;
}

struct Hit { int score; int64_t pos; char strand; std::string cigar; };

inline MapResult mapRead(const std::string& read, const KmerIndex& idx, const std::string& ref,
                         const MapParams& p, Stats& st) {
    using clk = std::chrono::steady_clock;
    auto secs = [](clk::time_point a, clk::time_point b) { return std::chrono::duration<double>(b - a).count(); };
    std::vector<Hit> hits; int ncand = 0;
    const std::string rc = revComp(read);
    for (int s = 0; s < 2; ++s) {                       // oba lanca: read i njegov reverse-complement
        const std::string& seq = (s == 0) ? read : rc;
        auto t0 = clk::now();
        auto cands = findCandidates(seq, idx, p);
        auto t1 = clk::now(); st.seedSec += secs(t0, t1);
        ncand += (int)cands.size();
        for (auto& c : cands) {
            int64_t rs = std::max<int64_t>(0, c.first - p.margin);
            int64_t re = std::min<int64_t>((int64_t)ref.size(), c.first + (int64_t)seq.size() + p.margin);
            if (re <= rs) continue;
            Alignment a = smithWaterman(seq, ref.substr((size_t)rs, (size_t)(re - rs)), p.sc);   // extension
            st.cells += (uint64_t)seq.size() * (uint64_t)(re - rs); ++st.extensions;
            if (a.score <= 0) continue;
            std::string cig;
            if (a.startX > 0) cig += std::to_string(a.startX) + "S";                         // soft clip
            cig += a.cigar;
            if ((int)seq.size() > a.endX) cig += std::to_string((int)seq.size() - a.endX) + "S";
            hits.push_back({a.score, rs + a.startY, s == 0 ? '+' : '-', cig});
        }
        st.extSec += secs(t1, clk::now());
    }
    MapResult r; r.candidates = ncand;
    if (hits.empty()) return r;
    std::sort(hits.begin(), hits.end(), [](const Hit& x, const Hit& y) { return x.score > y.score; });
    const Hit& b = hits[0];
    if (b.score < (int)(p.minScoreFrac * read.size() * p.sc.match)) return r;
    r.status = "mapped"; r.strand = b.strand; r.pos = b.pos; r.score = b.score; r.cigar = b.cigar;
    for (size_t i = 1; i < hits.size(); ++i)             // višeznačno: jednako dobar hit na drugom mjestu
        if (hits[i].score >= b.score && (hits[i].strand != b.strand || std::llabs(hits[i].pos - b.pos) > p.margin)) {
            r.status = "ambiguous"; break;
        }
    return r;
}

}  // namespace sw
