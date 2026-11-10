// gen_data: generator sintetičke reference i readova s poznatom pravom pozicijom.
// Upotreba: ./gen_data --ref-len 1000000 --reads 10000 --read-len 150 --sub 0.01 --ins 0.001 --del 0.001 --seed 42 --out data/sim
// Izlaz: <out>_ref.fa i <out>_reads.fq. ID reada: read_<i>|pos=<0-indeksirana pozicija>|strand=<+/->
#include <cstdio>
#include <iostream>
#include <map>
#include <random>
#include "fasta.hpp"
#include "mapper.hpp"
using namespace sw;

int main(int argc, char** argv) {
    std::map<std::string, std::string> o = {{"ref-len", "1000000"}, {"reads", "10000"}, {"read-len", "150"},
        {"sub", "0.01"}, {"ins", "0.001"}, {"del", "0.001"}, {"seed", "42"}, {"out", "data/sim"}};
    for (int i = 1; i + 1 < argc; i += 2) o[std::string(argv[i]).substr(2)] = argv[i + 1];
    const int L = std::stoi(o["ref-len"]), nReads = std::stoi(o["reads"]), rl = std::stoi(o["read-len"]);
    const double ps = std::stod(o["sub"]), pi = std::stod(o["ins"]), pd = std::stod(o["del"]);
    std::mt19937_64 rng(std::stoull(o["seed"]));
    std::uniform_real_distribution<double> U(0, 1);
    const char B[] = "ACGT";

    std::string ref(L, 'A');
    for (auto& c : ref) c = B[rng() & 3];
    std::vector<Record> reads;
    for (int i = 0; i < nReads; ++i) {
        int start = (int)(rng() % (uint64_t)(L - rl + 1));
        std::string r;
        for (char c : ref.substr(start, rl)) {
            if (U(rng) < pd) continue;                                     // delecija
            if (U(rng) < ps) { char n; do n = B[rng() & 3]; while (n == c); c = n; }   // zamjena
            r += c;
            if (U(rng) < pi) r += B[rng() & 3];                            // insercija
        }
        char strand = '+';
        if (rng() & 1) { r = revComp(r); strand = '-'; }
        char id[96]; std::snprintf(id, sizeof id, "read_%d|pos=%d|strand=%c", i, start, strand);
        reads.push_back({id, r, ""});
    }
    writeFasta(o["out"] + "_ref.fa", {{"ref", ref, ""}});
    writeFastq(o["out"] + "_reads.fq", reads);
    std::cerr << "Generisano: referenca " << L << " baza, " << nReads << " readova -> " << o["out"] << "_{ref.fa,reads.fq}\n";
}
