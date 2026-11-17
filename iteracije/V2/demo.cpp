// demo: mali primjer cijelog toka s ispisom ulaza i izlaza svakog koraka.
// Prevođenje: g++ -O2 -std=c++17 demo.cpp -o demo   (ili: make demo_small)
#include <iostream>
#include <random>
#include "mapper.hpp"
using namespace sw;
using std::cout;

static void showAlignment(const Alignment& a, const std::string& read, const std::string& region) {
    std::string mid;
    for (size_t i = 0; i < a.a.size(); ++i) mid += (a.a[i] == '-' || a.b[i] == '-') ? ' ' : (a.a[i] == a.b[i] ? '|' : 'x');
    cout << "    read : " << a.a << "   (iz reada koristeno [" << a.startX << "," << a.endX << ") od " << read.size() << ")\n"
         << "           " << mid << "\n"
         << "    ref  : " << a.b << "   (iz regije koristeno [" << a.startY << "," << a.endY << ") od " << region.size() << ")\n";
}

int main() {
    // ---- 1) ULAZ: referenca (nasumična, ali uvijek ista) ----
    std::mt19937 rng(7); const char B[] = "ACGT";
    std::string ref(300, 'A'); for (auto& c : ref) c = B[rng() & 3];
    cout << "=== ULAZ: referenca (" << ref.size() << " baza) ===\n" << ref.substr(0, 100) << "...\n\n";

    // ---- 2) ULAZ: readovi (izrezani iz reference, s ručno ubačenim greškama) ----
    struct R { std::string name, seq; std::string opis; };
    std::vector<R> reads;
    reads.push_back({"read_A", ref.substr(40, 50), "tačan isječak, pozicija 40"});
    { std::string s = ref.substr(120, 50); s[20] = (s[20] == 'A') ? 'C' : 'A';
      reads.push_back({"read_B", s, "pozicija 120, jedna zamjena na mjestu 20"}); }
    { std::string s = ref.substr(200, 50); s.erase(25, 1);
      reads.push_back({"read_C", s, "pozicija 200, jedna delecija na mjestu 25"}); }
    reads.push_back({"read_D", revComp(ref.substr(70, 50)), "pozicija 70, ali na suprotnom lancu (reverse-complement)"});
    { std::string s(50, 'A'); for (auto& c : s) c = B[rng() & 3];
      reads.push_back({"read_E", s, "nasumičan niz koji nije iz reference"}); }
    cout << "=== ULAZ: " << reads.size() << " readova ===\n";
    for (auto& r : reads) cout << "  " << r.name << "  " << r.seq << "\n      -> " << r.opis << "\n";

    // ---- 3) INDEKS ----
    MapParams p; p.k = 8; p.margin = 10; p.minVotes = 2;
    KmerIndex idx; idx.build(ref, p.k);
    cout << "\n=== INDEKS: k=" << p.k << ", " << idx.e.size() << " k-mera u indeksu ===\n";
    { std::string km = ref.substr(40, p.k); uint64_t code = 0; for (char c : km) code = (code << 2) | enc(c);
      auto r = idx.lookup(code);
      cout << "  lookup(" << km << ") -> pozicije:"; for (auto it = r.first; it != r.second; ++it) cout << ' ' << (uint32_t)(*it); cout << "\n"; }

    // ---- 4) MAPIRANJE svakog reada ----
    cout << "\n=== MAPIRANJE ===\n";
    Stats st; std::vector<MapResult> res;
    for (auto& r : reads) {
        cout << "\n" << r.name << " (" << r.opis << ")\n";
        for (int s = 0; s < 2; ++s) {
            std::string seq = s == 0 ? r.seq : revComp(r.seq);
            auto c = findCandidates(seq, idx, p);
            cout << "  lanac " << (s == 0 ? '+' : '-') << ": kandidati (dijagonala = pozicija početka, glasovi):";
            if (c.empty()) cout << " nema";
            for (auto& x : c) cout << "  (" << x.first << ", " << x.second << ")";
            cout << "\n";
            if (!c.empty()) {
                int64_t rs = std::max<int64_t>(0, c[0].first - p.margin), re = std::min<int64_t>(ref.size(), c[0].first + (int64_t)seq.size() + p.margin);
                std::string region = ref.substr(rs, re - rs);
                Alignment a = smithWaterman(seq, region, p.sc);
                cout << "    poravnanje najboljeg kandidata: skor " << a.score << ", CIGAR " << a.cigar << ", početak u referenci " << rs + a.startY << "\n";
                showAlignment(a, seq, region);
            }
        }
        res.push_back(mapRead(r.seq, idx, ref, p, st));
    }

    // ---- 5) IZLAZ: tabela rezultata (isti format kao results/map.tsv) ----
    cout << "\n=== IZLAZ ===\nread_id\tstrand\tpos\tscore\tcigar\tstatus\n";
    for (size_t i = 0; i < reads.size(); ++i)
        cout << reads[i].name << '\t' << res[i].strand << '\t' << res[i].pos << '\t' << res[i].score << '\t' << res[i].cigar << '\t' << res[i].status << '\n';
    cout << "\n(očekivano: A→40, B→120, C→200, D→70 na lancu -, E→unmapped)\n";
}
