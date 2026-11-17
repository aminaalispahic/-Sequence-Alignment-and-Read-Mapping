#!/usr/bin/env python3
"""Poredi Smith-Waterman skorove test_align (--stdin) s Biopython PairwiseAligner (local) na nasumičnim parovima.
Upotreba: pip install biopython && python3 scripts/compare_biopython.py [broj_parova]"""
import random, subprocess, sys
from Bio import Align

n = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
random.seed(7)
pairs = [("".join(random.choices("ACGT", k=random.randint(1, 80))),
          "".join(random.choices("ACGT", k=random.randint(1, 80)))) for _ in range(n)]

def aligner(mode):
    a = Align.PairwiseAligner()
    a.mode = mode
    a.match_score, a.mismatch_score = 3, -3
    a.open_gap_score = a.extend_gap_score = -2       # linearna kazna: 2 po znaku praznine
    return a

bad = 0
al = aligner("local")
inp = "".join(f"{x} {y}\n" for x, y in pairs)
out = subprocess.run(["./test_align", "--stdin"], input=inp, capture_output=True, text=True, check=True).stdout.split()
for (x, y), mine in zip(pairs, out):
    ref = int(al.score(x, y))
    if ref != int(mine):
        bad += 1; print(f"RAZLIKA: {x} {y} moj={mine} biopython={ref}")
print(f"Provjereno {n} poređenja, razlika: {bad}")
sys.exit(1 if bad else 0)
