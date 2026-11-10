#pragma once
#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

namespace sw {

struct Scoring { int match = 3, mismatch = -3, gap = 2; };

struct Alignment {
    int score = 0;
    std::string a, b;                                  
    int startX = 0, endX = 0, startY = 0, endY = 0;    
    std::string cigar;                                
};

enum Dir : uint8_t { STOP = 0, DIAG = 1, UP = 2, LEFT = 3 };

inline int sub(char x, char y, const Scoring& sc) {
    return (x == y && x != 'N') ? sc.match : sc.mismatch;  
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


inline Alignment smithWaterman(const std::string& X, const std::string& Y, const Scoring& sc) {
    const int n = (int)X.size(), m = (int)Y.size();
    const size_t W = (size_t)m + 1;
    std::vector<int> H((n + 1) * W, 0);                       
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


inline int smithWatermanScore(const std::string& X, const std::string& Y, const Scoring& sc) {
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


inline int rescore(const std::string& a, const std::string& b, const Scoring& sc) {
    int s = 0;
    for (size_t k = 0; k < a.size(); ++k)
        s += (a[k] == '-' || b[k] == '-') ? -sc.gap : sub(a[k], b[k], sc);
    return s;
}

}  // namespace sw
