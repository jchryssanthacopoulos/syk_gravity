#!/bin/bash
# BPS heterogeneity of the eta-averaged dressing cost (D7 sec. 7), N' = 10-16, two realizations each.
cd "$(dirname "$0")/../../.."
for N in 10 11 12 13 14 15 16; do
  .venv/bin/python scripts/analyze_bps_heterogeneity.py --Np $N --seeds 9000 9001
done
