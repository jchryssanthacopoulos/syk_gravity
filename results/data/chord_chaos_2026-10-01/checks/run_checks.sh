#!/bin/bash
# Reproduces the auxiliary checks of research/notes/chord_size_dictionary_2026-10-02.md (run from the repo root).
set -e
C=results/data/chord_chaos_2026-10-01/checks
PY=.venv/bin/python
$PY scripts/check_plethysm_rows.py --cases 2,3 3,3 2,5 --out $C/plethysm_rows.txt
$PY scripts/run_chord1_sizes.py --Np 10 --p 3 --samples 300 --seed 5 --out $C/chord1_p3_N10_F5.json
$PY scripts/run_chord1_sizes.py --Np 9 --F 4 --p 3 --samples 2000 --seed 5 --out $C/chord1_p3_N9_F4.json
$PY scripts/run_chord1_sizes.py --Np 10 --p 5 --samples 400 --seed 11 --out $C/chord1_p5_N10_F5.json
$PY scripts/run_chord1_sizes.py --Np 12 --p 5 --samples 300 --seed 11 --out $C/chord1_p5_N12_F6.json
rm -f $C/topsize_p3_N8.jsonl $C/topsize_p3_N10.jsonl $C/topsize_p3_N12.jsonl
$PY scripts/run_topsize_projection.py --Np 8 --nmax 2 --samples 2000 --seed 7 --out $C/topsize_p3_N8.jsonl | tail -1
$PY scripts/run_topsize_projection.py --Np 10 --nmax 2 --samples 1000 --seed 7 --out $C/topsize_p3_N10.jsonl | tail -1
$PY scripts/run_topsize_projection.py --Np 12 --nmax 3 --samples 300 --seed 21 --out $C/topsize_p3_N12.jsonl | tail -1
echo ALL_CHECKS_DONE
