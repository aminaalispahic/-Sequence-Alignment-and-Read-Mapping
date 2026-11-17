// mapper.hpp: k-mer indeks, seeding, extension (Smith-Waterman) i izbor najboljeg mapiranja
//
// Raspored podataka (savjeti s Laba 1):
//  - Readovi su u SoA obliku (ReadBatch): sve baze jedna iza druge u jednom nizu + niz početaka.
//    Umjesto vector<string> razbacanih po memoriji -> sekvencijalan pristup, dobra lokalnost.
//  - Baze se kodiraju jednom (0..3, N=4); referenca također.
//  - Indeks je direktna tabela (CSR): lookup = 2 čitanja iz memorije umjesto binarne pretrage.
#pragma once
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include "align.hpp"

namespace sw {

// Reverse-complement nad tekstom (koristi gen_data).
inline std::string revComp(const std::string& s) {
    std::string r(s.rbegin(), s.rend());
    for (char& c : r) switch (c) { case 'A': c = 'T'; break; case 'C': c = 'G'; break;
                                   case 'G': c = 'C'; break; case 'T': c = 'A'; break; default: break; }
    return r;
}
// Reverse-complement nad kodiranim bazama: komplement od b je 3-b (A<->T, C<->G), N ostaje N.
inline void revCompEnc(const uint8_t* s, int n, std::vector<uint8_t>& out) {
    out.resize(n);
    for (int i = 0; i < n; ++i) { const uint8_t b = s[n - 1 - i]; out[i] = (b == BASE_N) ? BASE_N : (uint8_t)(3 - b); }
}

// ---------------------------------------------------------------------------
// SoA skup readova: read i = bases[start[i] .. start[i+1]).
// ---------------------------------------------------------------------------
struct ReadBatch {
    std::vector<uint8_t> bases;
    std::vector<uint32_t> start{0};
    std::vector<std::string> ids;      // treba samo za ispis/evaluaciju, ne u vrućoj petlji
    size_t size() const { return ids.size(); }
    const uint8_t* data(size_t i) const { return bases.data() + start[i]; }
    int len(size_t i) const { return (int)(start[i + 1] - start[i]); }
};
template <class Records>
inline ReadBatch makeBatch(const Records& recs) {
    ReadBatch b;
    size_t total = 0; for (const auto& r : recs) total += r.seq.size();
    b.bases.reserve(total); b.start.reserve(recs.size() + 1); b.ids.reserve(recs.size());
    for (const auto& r : recs) {
        for (char c : r.seq) b.bases.push_back(enc(c));
        b.start.push_back((uint32_t)b.bases.size());
        b.ids.push_back(r.id);
    }
    return b;
}

// ---------------------------------------------------------------------------
// Indeks: direktna tabela (CSR). offs[code]..offs[code+1] su pozicije k-mera `code` u `pos`.
// Gradi se counting sortom: brojanje -> prefiksna suma -> popunjavanje (O(N), bez std::sort).
// Memorija: 4^k * 4 B (k=12 -> 64 MB), zato k <= 13.
// ---------------------------------------------------------------------------
struct KmerIndex {
    int k = 12;
    std::vector<uint32_t> offs, pos;
    void build(const std::vector<uint8_t>& ref, int k_) {
        k = k_;
        const uint64_t K = 1ULL << (2 * k), mask = K - 1;
        offs.assign(K + 1, 0);
        uint64_t code = 0; int valid = 0;
        for (size_t i = 0; i < ref.size(); ++i) {                  // 1) brojanje
            const uint8_t c = ref[i];
            if (c == BASE_N) { valid = 0; code = 0; continue; }
            code = ((code << 2) | c) & mask;
            if (++valid >= k) ++offs[code + 1];
        }
        for (uint64_t c = 0; c < K; ++c) offs[c + 1] += offs[c];  // 2) prefiksna suma
        pos.resize(offs[K]);
        std::vector<uint32_t> fill(offs.begin(), offs.end() - 1);
        code = 0; valid = 0;
        for (size_t i = 0; i < ref.size(); ++i) {                  // 3) popunjavanje (pozicije ostaju sortirane)
            const uint8_t c = ref[i];
            if (c == BASE_N) { valid = 0; code = 0; continue; }
            code = ((code << 2) | c) & mask;
            if (++valid >= k) pos[fill[code]++] = (uint32_t)(i - k + 1);
        }
    }
    std::pair<const uint32_t*, const uint32_t*> lookup(uint64_t code) const {
        return {pos.data() + offs[code], pos.data() + offs[code + 1]};
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
    double relVotes = 0.5;      // zadrži kandidate s >= relVotes * (najbolji broj glasova); 0 = isključeno
    bool scoreFirst = true;     // više kandidata: prvo SW samo skor, traceback samo za najbolje
    int band = 12;              // polu-širina trake za banded SW (-1 = cijeli prozor, tj. puni SW)
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
// k-mer se računa "rolling" (code zavisi od prethodnog) -> ta petlja se ne vektorizuje.
// Seeding je MEMORY BOUND (nasumični skokovi po indeksu od 64 MB). Zato tri prolaza sa software prefetchom:
// prvo se izračunaju svi k-meri i zatraže njihovi unosi, pa se tek onda čitaju. Procesor tako čeka na više
// keš promašaja istovremeno, a ne na jedan po jedan. Svi pomoćni nizovi se ponovo koriste (thread_local).
#if defined(__GNUC__) || defined(__clang__)
#define SW_PREFETCH(p) __builtin_prefetch(p)
#else
#define SW_PREFETCH(p) ((void)0)
#endif
inline const std::vector<std::pair<int64_t, int>>& findCandidates(const uint8_t* seq, int n, const KmerIndex& idx, const MapParams& p) {
    struct Cl { int64_t diag; int votes; int bestRun; };
    thread_local std::vector<uint32_t> codes, rs, re;
    thread_local std::vector<int> offs;
    thread_local std::vector<int64_t> diags;
    thread_local std::vector<Cl> cl;
    thread_local std::vector<std::pair<int64_t, int>> out;
    codes.clear(); offs.clear(); rs.clear(); re.clear(); diags.clear(); cl.clear(); out.clear();
    const int k = idx.k;
    const uint64_t mask = (1ULL << (2 * k)) - 1;
    uint64_t code = 0; int valid = 0;
    for (int i = 0; i < n; ++i) {                              // prolaz 1: k-meri + prefetch tabele
        const uint8_t c = seq[i];
        if (c == BASE_N) { valid = 0; code = 0; continue; }
        code = ((code << 2) | c) & mask;
        if (++valid < k) continue;
        const int off = i - k + 1;
        if (off % p.stride) continue;
        codes.push_back((uint32_t)code); offs.push_back(off);
        SW_PREFETCH(&idx.offs[code]);
    }
    for (uint32_t c : codes) {                                 // prolaz 2: opsezi + prefetch pozicija
        const uint32_t b = idx.offs[c], e = idx.offs[c + 1];
        rs.push_back(b); re.push_back(e);
        if (e > b) SW_PREFETCH(&idx.pos[b]);
    }
    for (size_t q = 0; q < codes.size(); ++q) {                // prolaz 3: dijagonale
        if ((int)(re[q] - rs[q]) > p.maxOcc) continue;
        for (uint32_t t = rs[q]; t < re[q]; ++t) diags.push_back((int64_t)idx.pos[t] - offs[q]);
    }
    std::sort(diags.begin(), diags.end());
    int64_t last = 0;
    for (size_t i = 0; i < diags.size();) {
        size_t j = i; while (j < diags.size() && diags[j] == diags[i]) ++j;
        const int cnt = (int)(j - i); const int64_t d = diags[i];
        if (cl.empty() || d - last > p.margin) cl.push_back({d, cnt, cnt});
        else { Cl& c = cl.back(); c.votes += cnt; if (cnt > c.bestRun) { c.bestRun = cnt; c.diag = d; } }
        last = d; i = j;
    }
    for (auto& c : cl) if (c.votes >= p.minVotes) out.push_back({c.diag, c.votes});
    std::sort(out.begin(), out.end(), [](auto& x, auto& y) { return x.second > y.second; });
    if ((int)out.size() > p.maxCand) out.resize(p.maxCand);
    return out;
}

struct Hit { int score; int64_t pos; char strand; std::string cigar; };

inline MapResult mapRead(const uint8_t* read, int n, const KmerIndex& idx, const std::vector<uint8_t>& ref,
                         const MapParams& p, Stats& st) {
    using clk = std::chrono::steady_clock;
    auto secs = [](clk::time_point a, clk::time_point b) { return std::chrono::duration<double>(b - a).count(); };
    thread_local std::vector<uint8_t> rc;
    revCompEnc(read, n, rc);
    const uint8_t* seqs[2] = {read, rc.data()};

    // --- seeding za oba lanca ---
    auto t0 = clk::now();
    struct Cand { int64_t diag; int votes; int strand; };
    std::vector<Cand> cands; int bestVotes = 0;
    for (int s = 0; s < 2; ++s)
        for (auto& c : findCandidates(seqs[s], n, idx, p)) { cands.push_back({c.first, c.second, s}); bestVotes = std::max(bestVotes, c.second); }
    // Nasumičan k-mer pogodak često da 2-3 glasa (susjedni k-meri se preklapaju), a pravi ~n-k.
    // Relativni prag izbaci lažne kandidate prije skupog poravnanja.
    if (p.relVotes > 0) {
        const int thr = std::max(p.minVotes, (int)(p.relVotes * bestVotes));
        cands.erase(std::remove_if(cands.begin(), cands.end(), [&](const Cand& c) { return c.votes < thr; }), cands.end());
    }
    auto t1 = clk::now(); st.seedSec += secs(t0, t1);

    // --- extension ---
    const int band = (p.band >= 0) ? p.band : (1 << 20);         // -1 = traka pokriva cijeli prozor
    struct Win { int score; int64_t rs, re, diag; int strand; };
    std::vector<Win> wins;
    const bool needScore = p.scoreFirst && cands.size() > 1;     // jedan kandidat: odmah traceback
    for (auto& c : cands) {
        const int64_t rs = std::max<int64_t>(0, c.diag - p.margin);
        const int64_t re = std::min<int64_t>((int64_t)ref.size(), c.diag + n + p.margin);
        if (re <= rs) continue;
        int score = -1;
        if (needScore) {
            const int m = (int)(re - rs);
            score = smithWatermanScoreBanded(seqs[c.strand], n, ref.data() + rs, m, (int)(c.diag - rs), band, p.sc);
            st.cells += (uint64_t)n * (uint64_t)std::min(m, 2 * band + 1);
            ++st.extensions;
        }
        wins.push_back({score, rs, re, c.diag, c.strand});
    }
    auto finish = [&](const Win& w) {                            // SW s tracebackom (CIGAR)
        const int m = (int)(w.re - w.rs);
        Alignment a = smithWatermanBanded(seqs[w.strand], n, ref.data() + w.rs, m, (int)(w.diag - w.rs), band, p.sc);
        st.cells += (uint64_t)n * (uint64_t)std::min(m, 2 * band + 1); ++st.extensions;
        std::string cig;
        if (a.startX > 0) cig += std::to_string(a.startX) + "S";            // soft clip
        cig += a.cigar;
        if (n > a.endX) cig += std::to_string(n - a.endX) + "S";
        return Hit{a.score, w.rs + a.startY, w.strand == 0 ? '+' : '-', cig};
    };
    std::vector<Hit> hits;
    if (!needScore) for (auto& w : wins) hits.push_back(finish(w));
    else {
        int best = 0; for (auto& w : wins) best = std::max(best, w.score);
        for (auto& w : wins) if (w.score == best && best > 0) hits.push_back(finish(w));   // traceback samo za najbolje
    }
    st.extSec += secs(t1, clk::now());

    MapResult r; r.candidates = (int)cands.size();
    hits.erase(std::remove_if(hits.begin(), hits.end(), [](const Hit& h) { return h.score <= 0; }), hits.end());
    if (hits.empty()) return r;
    std::sort(hits.begin(), hits.end(), [](const Hit& x, const Hit& y) { return x.score > y.score; });
    const Hit& b = hits[0];
    if (b.score < (int)(p.minScoreFrac * n * p.sc.match)) return r;
    r.status = "mapped"; r.strand = b.strand; r.pos = b.pos; r.score = b.score; r.cigar = b.cigar;
    for (size_t i = 1; i < hits.size(); ++i)             // višeznačno: jednako dobar hit na drugom mjestu
        if (hits[i].score >= b.score && (hits[i].strand != b.strand || std::llabs(hits[i].pos - b.pos) > p.margin)) {
            r.status = "ambiguous"; break;
        }
    return r;
}

}  // namespace sw