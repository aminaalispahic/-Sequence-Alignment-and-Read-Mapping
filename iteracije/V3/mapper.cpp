// mapper: sekvencijalni read mapper (k-mer indeks + seeding + Smith-Waterman extension)
// Upotreba:
//   ./mapper --ref data/sim_ref.fa --reads data/sim_reads.fq --out results/map.tsv
//            [--k 12] [--match 3] [--mismatch -3] [--gap 2] [--margin 20] [--min-votes 2]
//            [--max-cand 5] [--min-score 0.5] [--tol 10] [--eval 1]
// Izlaz (TSV): read_id  strand  pos  score  cigar  status   (pos je 0-indeksirana)
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <cstdint>
#include <cstdlib>
#include "fasta.hpp"
#include "mapper.hpp"

using namespace sw;
using clk = std::chrono::steady_clock;
static double since(clk::time_point t) { return std::chrono::duration<double>(clk::now() - t).count(); }

int main(int argc, char** argv) {
    std::map<std::string, std::string> o = {{"out", "results/map.tsv"}, {"k", "12"}, {"match", "3"}, {"mismatch", "-3"},
        {"gap", "2"}, {"margin", "20"}, {"min-votes", "2"}, {"max-cand", "5"}, {"min-score", "0.5"}, {"tol", "10"},
        {"max-occ", "500"}, {"eval", "1"}, {"profile", "0"}};
    for (int i = 1; i + 1 < argc; i += 2) o[std::string(argv[i]).substr(2)] = argv[i + 1];
    const bool profile = o["profile"] == "1";
    if (!o.count("ref") || !o.count("reads")) { std::cerr << "Potrebno: --ref <fasta> --reads <fastq/fasta>\n"; return 1; }

    MapParams p;
    p.k = std::stoi(o["k"]); p.margin = std::stoi(o["margin"]); p.minVotes = std::stoi(o["min-votes"]);
    p.maxCand = std::stoi(o["max-cand"]); p.minScoreFrac = std::stod(o["min-score"]); p.maxOcc = std::stoi(o["max-occ"]);
    p.sc = {std::stoi(o["match"]), std::stoi(o["mismatch"]), std::stoi(o["gap"])};
    if (p.k < 1 || p.k > 16) { std::cerr << "k mora biti u [1,16]\n"; return 1; }

    auto T0 = clk::now(); auto t = T0;
    std::string refBases;
    {
        const auto refs = readSequences(o["ref"]);
        if (refs.empty()) { std::cerr << "Prazna referenca\n"; return 1; }
        refBases = encodeSequence(refs[0].seq);
    }
    const EncodedView ref{refBases};
    const ReadBatch reads = makeBatch(readSequences(o["reads"]));
    double tLoad = since(t);

    t = clk::now(); KmerIndex idx; idx.build(ref, p.k); double tIndex = since(t);

    Stats st; std::vector<MapResult> res(reads.size());
    t = clk::now();
    for (size_t i = 0; i < reads.size(); ++i)
        res[i] = mapRead(reads.sequence(i), reads.sequence(i, true), idx, ref, p, st);
    double tMap = since(t);

    if (profile) {
        std::uint64_t checksum = 0;

        for (const auto& r : res) {
            checksum = checksum * 31
                + static_cast<std::uint64_t>(r.pos)
                + static_cast<std::uint64_t>(r.score);
        }

        std::fprintf(
            stderr,
            "Checksum: %llx\n",
            static_cast<unsigned long long>(checksum)
        );

        return 0;
    }

    t = clk::now();
    { std::ofstream f(o["out"]);
      if (!f) { std::cerr << "Ne mogu pisati u " << o["out"] << " (napravi direktorij)\n"; return 1; }
      f << "read_id\tstrand\tpos\tscore\tcigar\tstatus\n";
      for (size_t i = 0; i < reads.size(); ++i)
          f << reads.ids[i] << '\t' << res[i].strand << '\t' << res[i].pos << '\t' << res[i].score << '\t'
            << res[i].cigar << '\t' << res[i].status << '\n'; }
    double tOut = since(t), tTotal = since(T0);

    size_t nMapped = 0, nAmb = 0, nUn = 0, nCorrect = 0, nTruth = 0; const int tol = std::stoi(o["tol"]);
    for (size_t i = 0; i < reads.size(); ++i) {
        const auto& r = res[i];
        if (r.status == "unmapped") ++nUn; else if (r.status == "ambiguous") ++nAmb; else ++nMapped;
        const auto& id = reads.ids[i];
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
    std::fprintf(stderr, "Extension: %llu poravnanja, %.3e ćelija, %.4f GCUPS\n",
                 (unsigned long long)st.extensions, (double)st.cells, st.extSec > 0 ? st.cells / st.extSec / 1e9 : 0.0);
    return 0;
}
