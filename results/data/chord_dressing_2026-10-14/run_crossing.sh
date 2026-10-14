#!/bin/bash
# Crossed matter-chord diagram in the third moment (unbiased three-vector estimator), N' = 10-14.
cd "$(dirname "$0")/../../.."
.venv/bin/python scripts/run_chord_crossing.py --Np 10 --nsamp 2000 --seeds 9000 9001
.venv/bin/python scripts/run_chord_crossing.py --Np 11 --nsamp 1500 --seeds 9000 9001
.venv/bin/python scripts/run_chord_crossing.py --Np 12 --nsamp 1000 --seeds 9000 9001
.venv/bin/python scripts/run_chord_crossing.py --Np 13 --nsamp 600 --seeds 9000 9001
.venv/bin/python scripts/run_chord_crossing.py --Np 14 --nsamp 300 --seeds 9000 9001
