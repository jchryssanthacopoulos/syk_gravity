#!/bin/bash
# Second half of run_checks.sh (after the half-filling fix in scripts/run_chord1_sizes.py); run from the repo root.
set -e
C=results/data/chord_chaos_2026-10-01/checks
PY=.venv/bin/python
$PY scripts/run_chord1_sizes.py --Np 10 --p 5 --samples 400 --seed 11 --out $C/chord1_p5_N10_F5.json > /dev/null
$PY scripts/run_chord1_sizes.py --Np 12 --p 5 --samples 300 --seed 11 --out $C/chord1_p5_N12_F6.json > /dev/null
rm -f $C/topsize_p3_N8.jsonl $C/topsize_p3_N10.jsonl $C/topsize_p3_N12.jsonl
$PY scripts/run_topsize_projection.py --Np 8 --nmax 2 --samples 2000 --seed 7 --out $C/topsize_p3_N8.jsonl | tail -1
$PY scripts/run_topsize_projection.py --Np 10 --nmax 2 --samples 1000 --seed 7 --out $C/topsize_p3_N10.jsonl | tail -1
$PY scripts/run_topsize_projection.py --Np 12 --nmax 3 --samples 300 --seed 21 --out $C/topsize_p3_N12.jsonl | tail -1
echo ALL_CHECKS_DONE
