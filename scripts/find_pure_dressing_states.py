#!/usr/bin/env python3
r"""
Pure-dressing BPS states of the two-flavor pair uplift (docs/derivations.md D8 sec. 9.4).

Four-slot injectivity (a BPS state is determined by the old-harmonic parts of its four slots) holds for generic couplings
but fails on a codimension-2 set: where the reduced two-slot system M of the closed-form uplift is singular, a BPS state
with vanishing old-harmonic part in every slot exists.  This script finds such couplings by scaling two independent
families of the couplings that involve the new pair and minimizing sigma_min(M) in one sector, then checks the BPS
dimension (unchanged), the smallest bare fraction (-> 0) and the injectivity deficit (-> 1).

  generic model:      gE scales C[N,i,k], C[i,N,k] (k old) [the E couplings], gK scales C[i,j,N] [the K couplings]
  symmetrized model:  g1 scales C[N,i,k] (symmetric in i,k), g2 scales C[i,N,k] = C[i,k,N] (symmetry preserved)

Usage: .venv/bin/python scripts/find_pure_dressing_states.py --model generic --N 4 --seed 2 --sector 3 3
       .venv/bin/python scripts/find_pure_dressing_states.py --model sym --N 4 --seed 1 --sector 2 2
"""
import argparse
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from two_flavor import random_generic_couplings, random_symmetric_couplings  # noqa: E402
from two_flavor_uplift import PairUplift  # noqa: E402
from run_two_flavor_uplift import git_state  # noqa: E402

OUT = os.path.join(ROOT, "results/data/two_flavor_2026-10-23")


def scaled(C1, g, model):
    C = C1.copy()
    N = C.shape[0] - 1
    old = np.arange(N)
    if model == "generic":
        C[N][np.ix_(old, old)] *= g[0]
        C[np.ix_(old, [N], old)] *= g[0]
        C[np.ix_(old, old, [N])] *= g[1]
    else:
        C[N][np.ix_(old, old)] *= g[0]
        C[np.ix_(old, [N], old)] *= g[1]
        C[np.ix_(old, old, [N])] *= g[1]
        assert np.allclose(C, C.transpose(0, 2, 1))
    return C


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", choices=["sym", "generic"], required=True)
    ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--sector", type=int, nargs=2, required=True)
    ap.add_argument("--grid", type=int, default=19)
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    A, B = a.sector
    rng = np.random.default_rng(a.seed)
    C1 = (random_symmetric_couplings if a.model == "sym" else random_generic_couplings)(a.N + 1, rng)

    def smin(p):
        P = PairUplift(scaled(C1, np.exp(p), a.model))
        return np.linalg.svd(P.reduced_system(A, B)["M"], compute_uv=False).min(), P

    t0 = time.time()
    best = (np.inf, None)
    for l1 in np.linspace(-1.5, 3.0, a.grid):
        for l2 in np.linspace(-1.5, 3.0, a.grid):
            s, _ = smin((l1, l2))
            if s < best[0]:
                best = (s, (l1, l2))
    res = minimize(lambda p: np.log(smin(p)[0]), best[1], method="Nelder-Mead",
                   options=dict(xatol=1e-12, fatol=1e-12, maxiter=4000, maxfev=8000))
    s, P = smin(res.x)
    mu = P.bare_fraction_spectrum(A, B)
    rec = dict(model=a.model, N=a.N, seed=a.seed, sector=[A, B], scalings=[float(x) for x in np.exp(res.x)],
               grid_best_sigma_min=float(best[0]), sigma_min_M=float(s), d_new=P.d_new(A, B), hpl_dim=P.hpl_dim(A, B),
               d_new_unscaled=PairUplift(C1).d_new(A, B), min_bare_fraction=float(mu.min()),
               injectivity_deficit=int(P.injectivity_deficit(A, B)), git=git_state(),
               seconds=round(time.time() - t0, 1))
    if a.model == "sym" and A == B:
        rec["tower_bare_fraction"] = float(P.tower_slot_weights(A)[1])
    print(json.dumps(rec, indent=1))
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, f"pure_dressing_{a.model}_N{a.N}_s{a.seed}_{A}{B}.json"), "w") as fh:
        json.dump(rec, fh, indent=1)


if __name__ == "__main__":
    main()
