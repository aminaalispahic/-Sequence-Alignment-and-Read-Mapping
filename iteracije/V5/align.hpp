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
#include "bases.hpp"

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

inline int sub(Base x, Base y, const Scoring& sc) {
    return (x == y && x != BASE_N) ? sc.match : sc.mismatch;
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

struct SwWorkspace {
    std::vector<int> previous, current;
    std::vector<uint8_t> directions;
};

// Smith-Waterman: dva reda skorova i puna matrica smjerova za traceback.
template <class Sequence>
inline Alignment smithWatermanImpl(Sequence X, Sequence Y, const Scoring& sc,
                                   SwWorkspace& workspace) {
    const int n = (int)X.size(), m = (int)Y.size();
    const size_t width = (size_t)m + 1;

    auto& previous = workspace.previous;
    auto& current = workspace.current;
    auto& directions = workspace.directions;
    previous.resize(width);
    current.resize(width);
    directions.resize(((size_t)n + 1) * width);
    std::fill(previous.begin(), previous.end(), 0);

    int best = 0, bestI = 0, bestJ = 0;
    for (int i = 1; i <= n; ++i) {
        current[0] = 0;
        for (int j = 1; j <= m; ++j) {
            const int d = previous[j - 1] + sub(X[i - 1], Y[j - 1], sc);
            const int u = previous[j] - sc.gap;
            const int l = current[j - 1] - sc.gap;
            const int h = std::max(std::max(0, d), std::max(u, l));

            current[j] = h;
            directions[(size_t)i * width + j] =
                (h == 0) ? STOP : (h == d) ? DIAG : (h == u) ? UP : LEFT;
            if (h > best) {
                best = h;
                bestI = i;
                bestJ = j;
            }
        }
        previous.swap(current);
    }

    Alignment result;
    result.score = best;
    result.endX = bestI;
    result.endY = bestJ;
    int i = bestI, j = bestJ;
    while (i > 0 && j > 0) {
        const uint8_t direction = directions[(size_t)i * width + j];
        if (direction == STOP) break;

        if (direction == DIAG) {
            result.a += decodeBase(X[--i]);
            result.b += decodeBase(Y[--j]);
        } else if (direction == UP) {
            result.a += decodeBase(X[--i]);
            result.b += '-';
        } else {
            result.a += '-';
            result.b += decodeBase(Y[--j]);
        }
    }

    result.startX = i;
    result.startY = j;
    std::reverse(result.a.begin(), result.a.end());
    std::reverse(result.b.begin(), result.b.end());
    result.cigar = makeCigar(result.a, result.b);
    return result;
}

inline Alignment smithWaterman(EncodedView X, EncodedView Y, const Scoring& sc,
                               SwWorkspace& workspace) {
    return smithWatermanImpl(X, Y, sc, workspace);
}

inline Alignment smithWaterman(const std::string& X, std::string_view Y, const Scoring& sc,
                               SwWorkspace& workspace) {
    return smithWatermanImpl(std::string_view(X), Y, sc, workspace);
}

inline Alignment smithWaterman(const std::string& X, std::string_view Y, const Scoring& sc) {
    SwWorkspace workspace;
    return smithWaterman(X, Y, sc, workspace);
}

// Smith-Waterman samo skor: koristi iste redove, bez matrice smjerova.
template <class Sequence>
inline int smithWatermanScoreImpl(Sequence X, Sequence Y, const Scoring& sc,
                                  SwWorkspace& workspace) {
    const int n = (int)X.size(), m = (int)Y.size();
    const size_t width = (size_t)m + 1;
    auto& prev = workspace.previous;
    auto& cur = workspace.current;
    prev.resize(width);
    cur.resize(width);
    std::fill(prev.begin(), prev.end(), 0);

    int best = 0;
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

inline int smithWatermanScore(EncodedView X, EncodedView Y, const Scoring& sc,
                              SwWorkspace& workspace) {
    return smithWatermanScoreImpl(X, Y, sc, workspace);
}

inline int smithWatermanScore(const std::string& X, const std::string& Y, const Scoring& sc) {
    SwWorkspace workspace;
    return smithWatermanScoreImpl(std::string_view(X), std::string_view(Y), sc, workspace);
}

// Ponovno računanje skora iz poravnatih nizova (za samoprovjeru).
inline int rescore(const std::string& a, const std::string& b, const Scoring& sc) {
    int s = 0;
    for (size_t k = 0; k < a.size(); ++k)
        s += (a[k] == '-' || b[k] == '-') ? -sc.gap : sub(a[k], b[k], sc);
    return s;
}

}  // namespace sw
