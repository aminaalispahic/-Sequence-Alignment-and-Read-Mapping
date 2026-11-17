// mapper: sekvencijalni read mapper (k-mer indeks + seeding + Smith-Waterman extension)
// Upotreba:
//   ./mapper --ref data/sim_ref.fa --reads data/sim_reads.fq --out results/map.tsv
//            [--k 12] [--match 3] [--mismatch -3] [--gap 2] [--margin 20] [--min-votes 2]
//            [--max-cand 5] [--min-score 0.5] [--tol 10] [--eval 1]
//            [--band 12] [--rel-votes 0.5] [--score-first 1]
//            [--repeat 5] [--warmup 1]
// Izlaz (TSV): read_id  strand  pos  score  cigar  status   (pos je 0-indeksirana)
//
// Mjerenje (kao bench.hpp s Laba 1): faza mapiranja se pokreće `warmup` puta bez mjerenja,
// pa `repeat` puta s mjerenjem; ispisuje se medijan i minimum. Rezultat se uvijek "troši"
// (checksum + upis u fajl), pa kompajler ne može izbaciti posao kao mrtav kod.
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include "fasta.hpp"
#include "mapper.hpp"
using namespace sw;
using clk = std::chrono::steady_clock;
static double since(clk::time_point t) { return std::chrono::duration<double>(clk::now() - t).count(); }

int main(int argc, char** argv) {
    std::map<std::string, std::string> o = {{"out", "results/map.tsv"}, {"k", "12"}, {"match", "3"}, {"mismatch", "-3"},
        {"gap", "2"}, {"margin", "20"}, {"min-votes", "2"}, {"max-cand", "5"}, {"min-score", "0.5"}, {"tol", "10"},
        {"max-occ", "500"}, {"eval", "1"}, {"rel-votes", "0.5"}, {"score-first", "1"}, {"band", "12"},
        {"repeat", "1"}, {"warmup", "0"}};
    for (int i = 1; i + 1 < argc; i += 2) o[std::string(argv[i]).substr(2)] = argv[i + 1];
    if (!o.count("ref") || !o.count("reads")) { std::cerr << "Potrebno: --ref <fasta> --reads <fastq/fasta>\n"; return 1; }

    MapParams p;
    p.k = std::stoi(o["k"]); p.margin = std::stoi(o["margin"]); p.minVotes = std::stoi(o["min-votes"]);
    p.maxCand = std::stoi(o["max-cand"]); p.minScoreFrac = std::stod(o["min-score"]); p.maxOcc = std::stoi(o["max-occ"]);
    p.relVotes = std::stod(o["rel-votes"]); p.scoreFirst = o["score-first"] == "1"; p.band = std::stoi(o["band"]);
    p.sc = {std::stoi(o["match"]), std::stoi(o["mismatch"]), std::stoi(o["gap"])};
    const int repeat = std::max(1, std::stoi(o["repeat"])), warmup = std::max(0, std::stoi(o["warmup"]));
    if (p.k < 1 || p.k > 13) { std::cerr << "k mora biti u [1,13] (direktna tabela 4^k)\n"; return 1; }

    // --- učitavanje + kodiranje (jednom) ---
    auto T0 = clk::now(); auto t = T0;
    auto refs = readSequences(o["ref"]);
    if (refs.empty()) { std::cerr << "Prazna referenca\n"; return 1; }
    const std::vector<uint8_t> ref = encodeSeq(refs[0].seq);   // radi se s prvom sekvencom reference
    const ReadBatch reads = makeBatch(readSequences(o["reads"]));
    refs.clear();
    double tLoad = since(t);

    t = clk::now(); KmerIndex idx; idx.build(ref, p.k); double tIndex = since(t);

    // --- mapiranje: warmup + repeat mjerenja ---
    std::vector<MapResult> res(reads.size());
    Stats st; std::vector<double> runs; std::vector<Stats> runStats;
    uint64_t checksum = 0;
    for (int it = 0; it < warmup + repeat; ++it) {
        st = Stats{};
        t = clk::now();
        for (size_t i = 0; i < reads.size(); ++i) res[i] = mapRead(reads.data(i), reads.len(i), idx, ref, p, st);
        const double sec = since(t);
        checksum = 0; for (auto& r : res) checksum = checksum * 31 + (uint64_t)r.pos + (uint64_t)r.score;   // "consume"
        if (it >= warmup) { runs.push_back(sec); runStats.push_back(st); }
    }
    std::vector<double> sorted = runs; std::sort(sorted.begin(), sorted.end());
    const size_t mid = sorted.size() / 2;
    const double tMap = (sorted.size() % 2) ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2;
    // seeding/extension prikazujemo iz pokretanja čije je vrijeme najbliže medijanu (ne iz zadnjeg pokretanja)
    { size_t bi = 0; for (size_t q = 1; q < runs.size(); ++q) if (std::abs(runs[q] - tMap) < std::abs(runs[bi] - tMap)) bi = q;
      st = runStats[bi]; }

    t = clk::now();
    { std::ofstream f(o["out"]);
      if (!f) { std::cerr << "Ne mogu pisati u " << o["out"] << " (napravi direktorij)\n"; return 1; }
      f << "read_id\tstrand\tpos\tscore\tcigar\tstatus\n";
      for (size_t i = 0; i < reads.size(); ++i)
          f << reads.ids[i] << '\t' << res[i].strand << '\t' << res[i].pos << '\t' << res[i].score << '\t'
            << res[i].cigar << '\t' << res[i].status << '\n'; }
    double tOut = since(t);
    const double tTotal = tLoad + tIndex + tMap + tOut;   // jedno pokretanje (medijan mapiranja)

    size_t nMapped = 0, nAmb = 0, nUn = 0, nCorrect = 0, nTruth = 0; const int tol = std::stoi(o["tol"]);
    for (size_t i = 0; i < reads.size(); ++i) {
        const auto& r = res[i];
        if (r.status == "unmapped") ++nUn; else if (r.status == "ambiguous") ++nAmb; else ++nMapped;
        const std::string& id = reads.ids[i];
        size_t a = id.find("|pos="), b = id.find("|strand=");
        if (a != std::string::npos && b != std::string::npos) {
            ++nTruth;
            long tp = std::stol(id.substr(a + 5, b - a - 5)); char ts = id[b + 8];
            if (r.status != "unmapped" && r.strand == ts && std::llabs((long long)r.pos - tp) <= tol) ++nCorrect;
        }
    }
    std::fprintf(stderr, "Readova: %zu | mapped: %zu, ambiguous: %zu, unmapped: %zu\n", reads.size(), nMapped, nAmb, nUn);
    if (nTruth && o["eval"] == "1")
        std::fprintf(stderr, "Tačnost (poznata prava pozicija, tol=%d): %zu/%zu = %.2f%%\n", tol, nCorrect, nTruth, 100.0 * nCorrect / nTruth);
    std::fprintf(stderr, "Vrijeme [s]: učitavanje %.3f | indeks %.3f | mapiranje %.3f (seeding %.3f, extension %.3f) | ispis %.3f | ukupno %.3f\n",
                 tLoad, tIndex, tMap, st.seedSec, st.extSec, tOut, tTotal);
    if (repeat > 1)
        std::fprintf(stderr, "Mapiranje: %d mjerenja (+%d warmup), medijan %.3f s, min %.3f s, max %.3f s\n",
                     repeat, warmup, tMap, sorted.front(), sorted.back());
    std::fprintf(stderr, "Extension: %llu poravnanja, %.3e ćelija, %.4f GCUPS\n",
                 (unsigned long long)st.extensions, (double)st.cells, st.extSec > 0 ? st.cells / st.extSec / 1e9 : 0.0);
    std::fprintf(stderr, "Checksum: %llx\n", (unsigned long long)checksum);
    return 0;
}