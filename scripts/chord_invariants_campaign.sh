#!/usr/bin/env bash
# Exact runs for the 2026-09-28 chord-invariant campaign (action item 1).  Usage:
#   bash scripts/chord_invariants_campaign.sh A|B|C|E     (groups can run in parallel)
# Output: results/data/chord_invariants_2026-09-28/  (append-only JSONL, one line per realization)
set -euo pipefail
cd "$(dirname "$0")/.."
PY=.venv/bin/python
RUN="$PY scripts/run_chord_invariants.py"
OUT=results/data/chord_invariants_2026-09-28
# Every group from F on runs under the memory watchdog (CLAUDE.md): per-job cap, 30 GB project budget.
MW() { cap=$1; shift; $PY scripts/memwatch.py --limit-gb "$cap" --wait --log "$OUT/memwatch.jsonl" -- "$@"; }
case "${1:?group}" in
  A)  # q=3 small sizes, many realizations
      export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2
      for N in 5 6 7 8; do $RUN --q 3 --Np $N --seeds 0-59 --haar-seeds 0-29 --modes 5 --triples 3; done
      for N in 9 10;    do $RUN --q 3 --Np $N --seeds 0-39 --haar-seeds 0-19 --modes 6 --triples 4; done
      $RUN --q 3 --Np 11 --seeds 0-23 --haar-seeds 0-11 --modes 6 --triples 4 ;;
  B)  export OMP_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3 VECLIB_MAXIMUM_THREADS=3
      $RUN --q 3 --Np 12 --seeds 0-15 --haar-seeds 0-7 --modes 6 --triples 4
      $RUN --q 3 --Np 13 --seeds 0-9  --haar-seeds 0-4 --modes 6 --triples 4 ;;
  C)  export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 VECLIB_MAXIMUM_THREADS=4
      $RUN --q 3 --Np 14 --seeds 0-3 --haar-seeds 0-1 --modes 6 --triples 3
      $RUN --q 3 --Np 15 --seeds 0-2 --haar-seeds 0-1 --modes 5 --triples 2
      $RUN --q 3 --Np 16 --seeds 0-1 --haar-seeds 0   --modes 4 --triples 1 --method basis ;;
  E)  # q=5 (complement method for large N'); a ~ 0.87-0.99
      export OMP_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3 VECLIB_MAXIMUM_THREADS=3
      for N in 9 10 11 12 13 14; do $RUN --q 5 --Np $N --seeds 0-19 --haar-seeds 0-9 --modes 6 --triples 4; done
      $RUN --q 5 --Np 15 --seeds 0-7 --haar-seeds 0-3 --modes 6 --triples 3
      $RUN --q 5 --Np 16 --seeds 0-3 --haar-seeds 0-1 --modes 5 --triples 2
      $RUN --q 5 --Np 17 --seeds 0-1 --haar-seeds 0   --modes 4 --triples 1 ;;
  # Group E was stopped after q=5 N'=12 seed 14 (to parallelize); its remainder was run as E1, E2a, E2b:
  E1) export OMP_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3 VECLIB_MAXIMUM_THREADS=3
      $RUN --q 5 --Np 12 --seeds 15-19 --haar-seeds 0-9 --modes 6 --triples 4
      for N in 13 14; do $RUN --q 5 --Np $N --seeds 0-19 --haar-seeds 0-9 --modes 6 --triples 4; done ;;
  E2a) export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2
      $RUN --q 5 --Np 15 --seeds 0-7 --haar-seeds 0-3 --modes 6 --triples 3 ;;
  E2b) export OMP_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3 VECLIB_MAXIMUM_THREADS=3
      $RUN --q 5 --Np 16 --seeds 0-3 --haar-seeds 0-1 --modes 5 --triples 2
      $RUN --q 5 --Np 17 --seeds 0-1 --haar-seeds 0   --modes 4 --triples 1 ;;
  # The last lines of C (q=3 N'=16) and E2b (q=5 N'=17) were OOM-killed when run concurrently (old code: full
  # SVD, ~16 GB estimated, + an unbounded r x r cache, tens of GB).  Rerun with the fixed code (eigh null space,
  # bounded LRU cache; regression-identical on stored rows) and the same parameters, under memwatch:
  F1) export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 VECLIB_MAXIMUM_THREADS=4      # est. peak 9.1 GB
      MW 13 $RUN --q 3 --Np 16 --seeds 0-1 --haar-seeds 0 --modes 4 --triples 1 --method basis --cache-gb 3 ;;
  F2) export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 VECLIB_MAXIMUM_THREADS=4      # est. peak 7.1 GB
      MW 11 $RUN --q 5 --Np 17 --seeds 0-1 --haar-seeds 0 --modes 4 --triples 1 --cache-gb 4 ;;
esac
