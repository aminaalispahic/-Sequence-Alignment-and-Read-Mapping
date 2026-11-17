// align.hpp: sekvencijalni Smith-Waterman (lokalno poravnanje)
// Linearna kazna za prazninu: svaka praznina kosta `gap` (pozitivan broj).
// Konvencija: X = read (redovi matrice), Y = referenca (kolone), kao u SAM-u:
//   praznina u X -> 'D' (brisanje), praznina u Y -> 'I' (umetanje).
#pragma once
#include <algorithm>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace sw {

struct Scoring { int match = 3, mismatch = -3, gap = 2; };

struct Alignment {
    int score = 0;
    std::string a, b;                                  // poravnati nizovi ('-' = praznina)
    int startX = 0, endX = 0, startY = 0, endY = 0;    // 0-indeksirano, [start, end)
    std::string cigar;                                 // npr. 3M1D2M
};

enum Dir : uint8_t { STOP = 0, DIAG = 1, UP = 2, LEFT = 3 };

inline int sub(char x, char y, const Scoring& sc) {
    return (x == y && x != 'N') ? sc.match : sc.mismatch;   // N se tretira kao nepodudaranje
}

inline std::string makeCigar(const std::string& a, const std::string& b) {
    std::string out; char cur = 0; int cnt = 0;
    for (size_t k = 0; k < a.size(); ++k) {
        char op = (a[k] == '-') ? 'D' : (b[k] == '-') ? 'I' : 'M';
        if (op == cur) ++cnt;
        else { if (cnt) out += std::to_string(cnt) + cur; cur = op; cnt = 1; }
    }
    if (cnt) out += std::to_string(cnt) + cur;
    return out;
}

// Smith-Waterman: puna matrica + traceback od maksimuma do nule.
inline Alignment smithWaterman(const std::string& X, const std::string& Y, const Scoring& sc) {
    const int n = (int)X.size(), m = (int)Y.size();
    const size_t W = (size_t)m + 1;
    std::vector<int> H((n + 1) * W, 0);                       // prvi red i kolona = 0
    std::vector<uint8_t> T((n + 1) * W, STOP);
    int best = 0, bi = 0, bj = 0;
    for (int i = 1; i <= n; ++i)
        for (int j = 1; j <= m; ++j) {
            int d = H[(i-1)*W + j-1] + sub(X[i-1], Y[j-1], sc);
            int u = H[(i-1)*W + j] - sc.gap;
            int l = H[i*W + j-1] - sc.gap;
            int h = std::max(std::max(0, d), std::max(u, l));
            H[i*W + j] = h;
            T[i*W + j] = (h == 0) ? STOP : (h == d) ? DIAG : (h == u) ? UP : LEFT;
            if (h > best) { best = h; bi = i; bj = j; }
        }
    Alignment r; r.score = best; r.endX = bi; r.endY = bj;
    int i = bi, j = bj;
    while (i > 0 && j > 0 && H[i*W + j] > 0) {
        uint8_t t = T[i*W + j];
        if (t == DIAG)    { r.a += X[i-1]; r.b += Y[j-1]; --i; --j; }
        else if (t == UP) { r.a += X[i-1]; r.b += '-';    --i; }
        else              { r.a += '-';    r.b += Y[j-1]; --j; }
    }
    r.startX = i; r.startY = j;
    std::reverse(r.a.begin(), r.a.end()); std::reverse(r.b.begin(), r.b.end());
    r.cigar = makeCigar(r.a, r.b);
    return r;
}

// Smith-Waterman samo skor: dva reda, O(m) memorije.
inline int smithWatermanScore(std::string_view X, std::string_view Y, const Scoring& sc) {
    const int n = (int)X.size(), m = (int)Y.size();
    int best = 0;
    std::vector<int> prev(m + 1, 0), cur(m + 1, 0);
    for (int i = 1; i <= n; ++i) {
        cur[0] = 0;
        for (int j = 1; j <= m; ++j) {
            int h = std::max(std::max(0, prev[j-1] + sub(X[i-1], Y[j-1], sc)),
                             std::max(prev[j] - sc.gap, cur[j-1] - sc.gap));
            cur[j] = h; best = std::max(best, h);
        }
        std::swap(prev, cur);
    }
    return best;
}


// ---------------------------------------------------------------------------
// Kodirane baze: A=0, C=1, G=2, T=3, sve ostalo (N) = 4.
// Kodira se JEDNOM pri učitavanju; jezgra rade nad uint8_t nizovima (manje memorije,
// nema switch-a u petlji, a isti raspored kasnije koriste SIMD/CUDA verzije).
// ---------------------------------------------------------------------------
constexpr uint8_t BASE_N = 4;
inline uint8_t enc(char c) {
    switch (c) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; default: return BASE_N; }
}
inline char dec(uint8_t b) { return "ACGTN"[b > 4 ? 4 : b]; }
inline std::vector<uint8_t> encodeSeq(std::string_view s) {
    std::vector<uint8_t> v(s.size());
    for (size_t i = 0; i < s.size(); ++i) v[i] = enc(s[i]);
    return v;
}

// Banded Smith-Waterman, samo skor. Računa samo ćelije u traci širine 2*band+1 oko
// očekivane dijagonale j = i + diagOff. Dva reda H (stanu u L1 keš), baferi se ponovo koriste.
// Pristup je sekvencijalan (red po red, kolona po kolona) -> dobra prostorna lokalnost.
inline int smithWatermanScoreBanded(const uint8_t* __restrict__ X, int n, const uint8_t* __restrict__ Y, int m,
                                    int diagOff, int band, const Scoring& sc) {
    thread_local std::vector<int> bufA, bufB;
    if ((int)bufA.size() < m + 2) { bufA.resize(m + 2); bufB.resize(m + 2); }
    std::fill(bufA.begin(), bufA.begin() + m + 2, 0);
    std::fill(bufB.begin(), bufB.begin() + m + 2, 0);
    int* __restrict__ prev = bufA.data();
    int* __restrict__ cur = bufB.data();
    const int match = sc.match, mis = sc.mismatch, gap = sc.gap;
    int best = 0;
    for (int i = 1; i <= n; ++i) {
        const int jlo = std::max(1, i + diagOff - band), jhi = std::min(m, i + diagOff + band);
        if (jlo > jhi) { std::swap(prev, cur); continue; }
        cur[jlo - 1] = 0;                                  // lijevo od trake = 0
        const uint8_t xi = X[i - 1];
        // Ova petlja se NE može vektorizovati: cur[j] zavisi od cur[j-1] (kao prefix_scan).
        for (int j = jlo; j <= jhi; ++j) {
            const int s = (xi == Y[j - 1] && xi != BASE_N) ? match : mis;   // bez grananja (cmov)
            int h = prev[j - 1] + s;
            h = std::max(h, prev[j] - gap);
            h = std::max(h, cur[j - 1] - gap);
            h = std::max(h, 0);
            cur[j] = h;
            best = std::max(best, h);
        }
        std::swap(prev, cur);
    }
    return best;
}

// Banded Smith-Waterman s tracebackom. Smjer se pamti samo za ćelije u traci (n x (2*band+1) bajtova,
// row-major), a H samo u dva reda. Baferi su thread_local i ponovo se koriste.
inline Alignment smithWatermanBanded(const uint8_t* __restrict__ X, int n, const uint8_t* __restrict__ Y, int m,
                                     int diagOff, int band, const Scoring& sc) {
    band = std::min(band, n + m);                          // veća traka od matrice nema smisla
    const int W = 2 * band + 1;
    thread_local std::vector<int> bufA, bufB;
    thread_local std::vector<uint8_t> T;
    if ((int)bufA.size() < m + 2) { bufA.resize(m + 2); bufB.resize(m + 2); }
    if (T.size() < (size_t)(n + 1) * W) T.resize((size_t)(n + 1) * W);
    std::fill(bufA.begin(), bufA.begin() + m + 2, 0);
    std::fill(bufB.begin(), bufB.begin() + m + 2, 0);
    int* prev = bufA.data(); int* cur = bufB.data();
    int best = 0, bi = 0, bj = 0;
    for (int i = 1; i <= n; ++i) {
        const int base = i + diagOff - band;               // kolona koja odgovara T[i][0]
        const int jlo = std::max(1, base), jhi = std::min(m, i + diagOff + band);
        uint8_t* Ti = T.data() + (size_t)i * W;
        std::fill(Ti, Ti + W, (uint8_t)STOP);
        if (jlo > jhi) { std::swap(prev, cur); continue; }
        cur[jlo - 1] = 0;
        const uint8_t xi = X[i - 1];
        for (int j = jlo; j <= jhi; ++j) {
            const int d = prev[j - 1] + ((xi == Y[j - 1] && xi != BASE_N) ? sc.match : sc.mismatch);
            const int u = prev[j] - sc.gap, l = cur[j - 1] - sc.gap;
            const int h = std::max(std::max(0, d), std::max(u, l));
            cur[j] = h;
            Ti[j - base] = (h == 0) ? STOP : (h == d) ? DIAG : (h == u) ? UP : LEFT;
            if (h > best) { best = h; bi = i; bj = j; }
        }
        std::swap(prev, cur);
    }
    Alignment r; r.score = best; r.endX = bi; r.endY = bj;
    int i = bi, j = bj;
    while (i > 0 && j > 0) {
        const int idx = j - (i + diagOff - band);
        if (idx < 0 || idx >= W) break;
        const uint8_t t = T[(size_t)i * W + idx];
        if (t == STOP) break;
        if (t == DIAG)    { r.a += dec(X[i-1]); r.b += dec(Y[j-1]); --i; --j; }
        else if (t == UP) { r.a += dec(X[i-1]); r.b += '-';         --i; }
        else              { r.a += '-';         r.b += dec(Y[j-1]); --j; }
    }
    r.startX = i; r.startY = j;
    std::reverse(r.a.begin(), r.a.end()); std::reverse(r.b.begin(), r.b.end());
    r.cigar = makeCigar(r.a, r.b);
    return r;
}

// Ponovno računanje skora iz poravnatih nizova (za samoprovjeru).
inline int rescore(const std::string& a, const std::string& b, const Scoring& sc) {
    int s = 0;
    for (size_t k = 0; k < a.size(); ++k)
        s += (a[k] == '-' || b[k] == '-') ? -sc.gap : sub(a[k], b[k], sc);
    return s;
}

}  // namespace sw