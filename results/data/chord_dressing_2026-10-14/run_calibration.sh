#!/bin/bash
# Calibration data for the chord dressing analysis (N' = 10-16), regenerated with the committed scripts.
cd "$(dirname "$0")/../../.."
for N in 10 11 12 13 14; do
  .venv/bin/python scripts/run_chord_dressing.py --Np $N --seeds 9000 9001 9002 9003 --neta 24
done
.venv/bin/python scripts/run_chord_dressing.py --Np 15 --seeds 9000 9001 9002 --neta 12
.venv/bin/python scripts/run_chord_dressing.py --Np 16 --seeds 9000 9001 --neta 6
