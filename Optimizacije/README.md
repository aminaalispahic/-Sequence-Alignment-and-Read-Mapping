# Sekvencijalni Smith–Waterman i read mapper (P21)

## Prevođenje i pokretanje
Kompajler: **g++ (GCC) ili clang++**, ne MSVC. Na Windowsu: MSYS2 (`mingw32-make` umjesto `make`).

```
make            # release build (-O2)
make test       # testovi
make demo       # mali primjer
make bench      # 5 mil. baza, 50.000 readova, 5 mjerenja + 1 warmup
```

Bez `make`:
```
g++ -O2 -std=c++17 -fno-tree-vectorize -o mapper mapper.cpp
g++ -O2 -std=c++17 -o gen_data gen_data.cpp
g++ -O2 -std=c++17 -fno-tree-vectorize -o test_align test_align.cpp
```

| Fajl | Sadržaj |
|---|---|
| `align.hpp` | SW (puni + traceback + CIGAR), SW score-only, **banded SW** (score-only i s tracebackom) nad kodiranim bazama |
| `fasta.hpp` | Čitanje FASTA/FASTQ (auto-detekcija), pisanje FASTA/FASTQ |
| `mapper.hpp` | SoA readovi (`ReadBatch`), k-mer indeks (direktna tabela), seeding s glasanjem po dijagonali, extension, izbor najboljeg hita |
| `mapper.cpp` | CLI mapper, mjerenje (warmup + medijan), tačnost, GCUPS |
| `gen_data.cpp` | Generator sintetičke reference i readova |
| `test_align.cpp` | Jedinični i nasumični testovi (uklj. banded == puni SW); `--stdin` za poređenje s Biopythonom |
| `compare_biopython.py` | Poređenje skorova s Biopython `PairwiseAligner` |

## Optimizacije serijske verzije (prema Labu 1)
| Optimizacija | Zašto |
|---|---|
| Release build `-O2`, GCC/Clang | Optimizacije kompajlera; MSVC slabo vektorizuje |
| `-fno-tree-vectorize` za serijsku verziju | Pošten baseline za poređenje sa SIMD/OpenMP/CUDA |
| SoA readovi: jedan niz baza + niz početaka | Sekvencijalan pristup memoriji umjesto `vector<string>` |
| Baze kodirane jednom (A=0,C=1,G=2,T=3,N=4) | Nema ponovnog kodiranja u petljama, manje memorije |
| Indeks kao direktna tabela (counting sort + prefiksna suma) | Seeding je memory bound: 2 čitanja umjesto binarne pretrage (~23 keš promašaja) |
| Relativni prag glasova (`--rel-votes 0.5`) | Izbacuje lažne kandidate (~10 → ~1 poravnanje po readu) |
| Banded SW (`--band 12`), H u 2 reda | Manje ćelija; radni skup stane u L1 keš → compute bound |
| Traceback samo za najbolji kandidat | Ostali kandidati: samo skor |
| Seeding u 3 prolaza sa software prefetchom | Seeding je memory bound: procesor čeka više keš promašaja istovremeno |
| `thread_local` baferi, bez `substr` kopija | Nema alokacija po readu (spremno za OpenMP) |
| Warmup + medijan, checksum rezultata | Pouzdano mjerenje; kompajler ne može izbaciti posao kao mrtav kod |

Zašto se SW petlja ne vektorizuje: `cur[j]` zavisi od `cur[j-1]` (kao `prefix_scan`).
Kompajler javlja `possible dependence` (`-fopt-info-vec-missed`). Rješenje u paralelnoj fazi:
SIMD preko više readova (SoA), anti-dijagonale, tiling za duge nizove.

## Ulaz / izlaz mappera
- Ulaz: `--ref ref.fa` (FASTA), `--reads reads.fq` (FASTQ ili FASTA).
- Parametri: `--k 12 --match 3 --mismatch -3 --gap 2 --margin 20 --min-votes 2 --max-cand 5 --min-score 0.5 --tol 10`
  `--band 12 --rel-votes 0.5 --score-first 1 --repeat 1 --warmup 0`.
- Stari način (za poređenje u izvještaju): `--band -1 --rel-votes 0 --score-first 0`.
- Izlaz: TSV `read_id strand pos score cigar status` (pos 0-indeksirana; status: mapped / ambiguous / unmapped).
- Na stderr: broj mapiranih, tačnost (ako ID reada sadrži `|pos=..|strand=..`), vremena po fazama, GCUPS, checksum.
- k mora biti ≤ 13 (tabela ima 4^k ulaza; k=12 → 64 MB).

## Konvencije
- X = read (redovi), Y = referenca (kolone); CIGAR kao u SAM-u (`D` = praznina u readu, `I` = praznina u referenci, `S` = soft clip).
- Kazna za prazninu `--gap` je pozitivan broj (linearna), `--mismatch` je negativan.
- `N` se tretira kao nepodudaranje, k-mere s `N` se ne indeksiraju.