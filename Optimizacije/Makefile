# Release build (optimizacije uključene). Kompajler: g++ (GCC) ili clang++, ne MSVC.
#   make            -> gen_data, mapper, test_align
#   make test       -> pokreće testove
#   make demo       -> mali primjer (generisanje + mapiranje)
#   make bench      -> veći primjer, 5 mjerenja + 1 warmup (medijan)
# Serijska verzija se prevodi s -fno-tree-vectorize (kao skalarna verzija na Labu 1),
# da bi poređenje sa SIMD/OpenMP/CUDA verzijama kasnije bilo pošteno.
CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra
SERIAL   := -fno-tree-vectorize

ifeq ($(OS),Windows_NT)
  EXE := .exe
endif

all: gen_data$(EXE) mapper$(EXE) test_align$(EXE)

mapper$(EXE): mapper.cpp mapper.hpp align.hpp fasta.hpp
	$(CXX) $(CXXFLAGS) $(SERIAL) -o $@ mapper.cpp
gen_data$(EXE): gen_data.cpp mapper.hpp align.hpp fasta.hpp
	$(CXX) $(CXXFLAGS) -o $@ gen_data.cpp
test_align$(EXE): test_align.cpp align.hpp
	$(CXX) $(CXXFLAGS) $(SERIAL) -o $@ test_align.cpp

test: test_align$(EXE)
	./test_align$(EXE)

dirs:
	@mkdir -p data results

demo: all dirs
	./gen_data$(EXE) --ref-len 1000000 --reads 10000 --out data/sim
	./mapper$(EXE) --ref data/sim_ref.fa --reads data/sim_reads.fq --out results/sim.tsv

bench: all dirs
	./gen_data$(EXE) --ref-len 5000000 --reads 50000 --sub 0.02 --out data/big
	./mapper$(EXE) --ref data/big_ref.fa --reads data/big_reads.fq --out results/big.tsv --repeat 5 --warmup 1

clean:
	rm -f gen_data$(EXE) mapper$(EXE) test_align$(EXE)

.PHONY: all test dirs demo bench clean