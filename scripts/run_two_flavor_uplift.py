#!/usr/bin/env python3
r"""
Exact pair uplift N -> N+1 of the two-flavor N=2 SYK model (docs/derivations.md D8; src/two_flavor_uplift.py).

Tasks (any subset, one JSONL record per realization and sector):
  dims    two-step dimension formula (Theorem A twice) and four-slot injectivity             [per new sector]
  cv      CV channel decoders: rank, constrained room, exact-lift locus vs generic count, U_T (CV eq. 89)
                                                                                              [per old sector]
  bare    bare states (kernel of the bare map) vs the generic count; bare-fraction spectrum quantiles;
          tower slot weights                                                                  [per new sector]
  first   exact vs first-order bare-fraction spectrum; incoherent vs interference first-order dressing
                                                                                              [per new sector]
  tower   CV dressing costs kappa_00, kappa_11 of the tower vs CV's bounds (eqs. 87-88)       [per rung n]

Couplings: rng = default_rng(seed); C = random_symmetric_couplings(N+1, rng) (or random_generic_couplings with
--model generic).  The old model is the restriction of the same couplings to the first N pairs.

Usage:  .venv/bin/python scripts/run_two_flavor_uplift.py --N 2 3 4 5 6 --seeds 1 2 3 4 5 --tasks cv
        .venv/bin/python scripts/run_two_flavor_uplift.py --N 7 --dry-run
Run N >= 6 under scripts/memwatch.py.
"""
import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from two_flavor import random_generic_couplings, random_symmetric_couplings  # noqa: E402
from two_flavor_uplift import CHANNELS, PairUplift, memory_estimate_gb  # noqa: E402

OUT = os.path.join(ROOT, "results/data/two_flavor_2026-10-23")


def git_state():
    try:
        g = subprocess.run(["git", "-C", ROOT, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", ROOT, "status", "--porcelain", "src", "scripts"], capture_output=True,
                               text=True).stdout.strip()
        return g + ("-dirty" if dirty else "")
    except OSError:
        return "unknown"


def q(x, ps=(0.1, 0.25, 0.5, 0.75, 0.9)):
    return [float(v) for v in np.quantile(x, ps)] if len(x) else []


def task_dims(P):
    out = []
    for A in range(P.N + 2):
        for B in range(P.N + 2):
            if P.new_dim(A, B) == 0:
                continue
            r = P.two_step_dimensions(A, B)
            r.update(A=A, B=B, injectivity_deficit=P.injectivity_deficit(A, B))
            out.append(r)
    return out


def task_cv(P, symmetrized):
    out = []
    for a in range(P.N + 1):
        for b in range(P.N + 1):
            r = P.cv_channels(a, b, symmetrized=symmetrized)
            if r is not None:
                out.append(r)
    return out


def task_bare(P, symmetrized):
    out = []
    for A in range(P.N + 2):
        for B in range(P.N + 2):
            d = P.d_new(A, B)
            if d == 0:
                continue
            g = P.bare_count_generic(A, B)
            gx = P.bare_count_generic(A, B, exact_subkernels=True)
            mu = P.bare_fraction_spectrum(A, B)
            r = dict(A=A, B=B, d_new=d, bare=P.bare_count(A, B), gen=g["gen"], gen_exact_sub=gx["gen"],
                     n01=g["n01"], n10=g["n10"], n01_exact=gx["n01"], n10_exact=gx["n10"], slots=g["slots"],
                     rooms=g["rooms"], n_mu_one=int(np.sum(mu > 1 - 1e-9)), mu_min=float(mu.min()),
                     mu_mean=float(mu.mean()), mu_q=q(mu))
            if symmetrized and A == B:
                w, bare = P.tower_slot_weights(A)
                r.update(tower_slot_weights=w, tower_bare_fraction=bare)
            out.append(r)
    return out


def task_first(P, symmetrized):
    out = []
    for A in range(P.N + 2):
        for B in range(P.N + 2):
            r = P.first_order(A, B)
            if r is None:
                continue
            me, m1, k1, ki = r["mu_exact"], r["mu_first"], r["kappa1"], r["kappa_inc"]
            dressed = me < 1 - 1e-8
            frac = (k1[dressed] - ki[dressed]) / ki[dressed] if dressed.any() else np.zeros(0)
            rec = dict(A=A, B=B, d_new=len(me), mu_exact_mean=float(me.mean()), mu_first_mean=float(m1.mean()),
                       mu_exact_q=q(me), mu_first_q=q(m1), n_one_exact=int(np.sum(me > 1 - 1e-9)),
                       n_one_first=int(np.sum(m1 > 1 - 1e-9)),
                       H_exact=float(np.mean((2 * me - 1) ** 2 * (1 - (2 * me - 1) ** 2))),
                       H_first=float(np.mean((2 * m1 - 1) ** 2 * (1 - (2 * m1 - 1) ** 2))),
                       interference_over_incoherent_q=q(frac), interference_over_incoherent_mean=float(frac.mean()) if len(frac) else None,
                       mu_exact=me.tolist(), mu_first=m1.tolist(), interference_over_incoherent=frac.tolist())
            out.append(rec)
    return out


def task_tower(P):
    out = []
    for n in range(0, P.N + 1):
        for ch in ("00", "11"):
            if ch == "00" and n == 0:
                continue
            out.append(P.tower_cost(n, ch))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--N", type=int, nargs="+", required=True, help="old size(s) N (uplift N -> N+1)")
    ap.add_argument("--seeds", type=int, nargs="+", default=[1])
    ap.add_argument("--model", choices=["sym", "generic"], default="sym")
    ap.add_argument("--tasks", nargs="+", default=["dims", "cv", "bare", "first", "tower"],
                    choices=["dims", "cv", "bare", "first", "tower"])
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--tag", default="")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    for N in a.N:
        est = memory_estimate_gb(N)
        print(f"N={N}->{N + 1}: largest new sector {est['dmax']}, estimated peak {est['total']:.2f} GB", flush=True)
    if a.dry_run:
        return
    os.makedirs(a.out, exist_ok=True)
    stamp = time.strftime("%Y%m%dT%H%M%S")
    meta = dict(args=vars(a), git=git_state(), started=stamp,
                memory_estimate_gb={N: memory_estimate_gb(N)["total"] for N in a.N})
    with open(os.path.join(a.out, f"meta_{a.model}_{'-'.join(map(str, a.N))}_{'-'.join(a.tasks)}{a.tag}_{stamp}.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    sym = a.model == "sym"
    for N in a.N:
        for seed in a.seeds:
            t0 = time.time()
            rng = np.random.default_rng(seed)
            P = PairUplift((random_symmetric_couplings if sym else random_generic_couplings)(N + 1, rng))
            for task in a.tasks:
                t1 = time.time()
                recs = dict(dims=lambda: task_dims(P), cv=lambda: task_cv(P, sym), bare=lambda: task_bare(P, sym),
                            first=lambda: task_first(P, sym), tower=lambda: task_tower(P))[task]()
                path = os.path.join(a.out, f"{task}_N{N}_{a.model}{a.tag}.jsonl")
                with open(path, "a") as fh:
                    for r in recs:
                        fh.write(json.dumps(dict(N=N, seed=seed, model=a.model, **r)) + "\n")
                print(f"N={N} seed={seed} {task}: {len(recs)} records [{time.time() - t1:.0f}s]", flush=True)
            print(f"N={N} seed={seed} done [{time.time() - t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
