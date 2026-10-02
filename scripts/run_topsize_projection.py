#!/usr/bin/env python3
r"""
Top-size (diagonal) weights of the chord sectors of the BPS projector, p = 3, half filling.

For n >= 1 the size-2n part of H^n is pure Wiener chaos n (a chaos-m component carries size <= 2m), so
X_n = (H^n)^(2n) is an exact chaos-n, size-2n direction and, by Cauchy-Schwarz,
        J_{n,2n} >= J^(H)_{n,2n} = (E <P, X_n>)^2 / (d E ||X_n||^2).
At N' = 8-12 this bound is saturated to <= 1% (n = 2) and ~5% (n = 3) by the OU fits of scripts/run_chord_chaos.py
(research/notes/chord_size_dictionary_2026-10-02.md, Sec. 6).  Per-realization cost: two dense products and n
Casimir applications per power (src/chord_chaos.TopSizeFilter); memory O(1) dense D x D matrices.

Usage (N' = 14: ~20 s per realization unloaded, peak ~0.6 GB):
  .venv/bin/python scripts/memwatch.py --limit-gb 3 --log results/data/chord_chaos_2026-10-01/topsize/memwatch.jsonl -- \
      .venv/bin/python scripts/run_topsize_projection.py --Np 14 --samples 40 --seed 41 \
      --out results/data/chord_chaos_2026-10-01/topsize/topsize_p3_N14.jsonl
Output: one json line per realization (overlaps <P,X_n> and norms ||X_n||^2); running estimates on stdout.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from chord_chaos import Sector, TopSizeFilter, bps_projector_hodge, random_couplings  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Np", type=int, required=True, help="even N' (half filling)")
    ap.add_argument("--nmax", type=int, default=3, help="largest chord sector n (needs 2n <= N'/2)")
    ap.add_argument("--samples", type=int, default=40)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    S = Sector(a.Np, 3)
    filt = TopSizeFilter(a.Np, S.F)
    ns = [n for n in range(2, a.nmax + 1) if 2 * n <= S.F]
    rng = np.random.default_rng(a.seed)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    acc = []
    t0 = time.time()
    for s in range(a.samples):
        c = random_couplings(rng, S.M)
        H = S.hamiltonian(c)
        P, _, d = bps_projector_hodge(S, c)
        Hn = H
        r = []
        for n in range(2, max(ns) + 1):
            Hn = Hn @ H
            if n in ns:
                X = filt.component(Hn, 2 * n)
                r += [float(np.vdot(P, X).real), float(np.vdot(X, X).real)]
                del X
        del H, P, Hn
        acc.append(r)
        with open(a.out, "a") as fh:
            fh.write(json.dumps(dict(Np=a.Np, seed=a.seed, s=s, d=int(d), ns=ns, r=r)) + "\n")
        m = np.array(acc).mean(0)
        est = {f"J{n}{2*n}": m[2 * i] ** 2 / (d * m[2 * i + 1]) for i, n in enumerate(ns)}
        print(json.dumps(dict(n_real=len(acc), secs=round(time.time() - t0, 1), **est)), flush=True)


if __name__ == "__main__":
    main()
