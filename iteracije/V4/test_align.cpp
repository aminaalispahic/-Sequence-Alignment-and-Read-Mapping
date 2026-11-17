// test_align: testovi korektnosti Smith-Waterman algoritma.
//   ./test_align            pokreće testove
//   ./test_align --stdin    čita linije "<X> <Y>" i ispisuje SW skor (za poređenje s Biopythonom)
#include <cstdio>
#include <iostream>
#include <random>
#include "align.hpp"
using namespace sw;

static int fails = 0;
#define CHECK(cond, msg) do { if (!(cond)) { ++fails; std::printf("FAIL: %s\n", msg); } } while (0)

int main(int argc, char** argv) {
    Scoring sc;                                       // +3 / -3 / gap 2
    if (argc > 1 && std::string(argv[1]) == "--stdin") {
        std::string x, y;
        while (std::cin >> x >> y) std::printf("%d\n", smithWatermanScore(x, y, sc));
        return 0;
    }
    // 1) demo iz plana
    auto r = smithWaterman("TGTTACGG", "GGTTGACTA", sc);
    CHECK(r.score == 13, "SW demo skor 13");
    CHECK(r.a == "GTT-AC" && r.b == "GTTGAC", "SW demo poravnanje");
    CHECK(r.cigar == "3M1D2M", "SW demo CIGAR");
    CHECK(r.startX == 1 && r.endX == 6 && r.startY == 1 && r.endY == 7, "SW demo pozicije");
    // 2) primjer iz teksta: AGC / AC, +3 / -2 / gap 1 -> 5
    Scoring ex{3, -2, 1};
    CHECK(smithWaterman("AGC", "AC", ex).score == 5, "SW AGC/AC skor 5");
    CHECK(smithWaterman("AGC", "AC", ex).a == "AGC" && smithWaterman("AGC", "AC", ex).b == "A-C", "SW AGC/AC poravnanje");
    // 3) identični nizovi, potpuno različiti nizovi, prazni nizovi
    CHECK(smithWaterman("ACGTACGT", "ACGTACGT", sc).score == 24, "SW identični");
    CHECK(smithWaterman("AAAA", "CCCC", sc).score == 0, "SW bez sličnosti = 0");
    CHECK(smithWaterman("", "ACG", sc).score == 0 && smithWaterman("ACG", "", sc).score == 0, "SW prazan niz");
    // 4) nasumični testovi: puna i score-only verzija, ponovni obračun skora, pozicije poravnanja
    std::mt19937 rng(1); const char B[] = "ACGT";
    for (int it = 0; it < 3000; ++it) {
        std::string x(1 + rng() % 60, 'A'), y(1 + rng() % 60, 'A');
        for (auto& c : x) c = B[rng() & 3];
        for (auto& c : y) c = B[rng() & 3];
        auto l = smithWaterman(x, y, sc);
        CHECK(l.score == smithWatermanScore(x, y, sc), "SW full == score-only");
        CHECK(rescore(l.a, l.b, sc) == l.score, "SW rescore");
        CHECK(l.score >= 0, "SW >= 0");
        std::string ra, rb;
        for (char c : l.a) if (c != '-') ra += c;
        for (char c : l.b) if (c != '-') rb += c;
        CHECK(ra == x.substr(l.startX, l.endX - l.startX) && rb == y.substr(l.startY, l.endY - l.startY), "SW pozicije");
    }
    std::printf(fails ? "NEUSPJEH: %d grešaka\n" : "Svi testovi uspješni.\n", fails);
    return fails ? 1 : 0;
}
