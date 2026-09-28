#!/usr/bin/env python3
r"""
Measure chord invariants of BPS-projected slot projectors (action item 1).  See src/chord_invariants.py.

For each realization (SYK couplings seed s, or a Haar projector of the same rank) this records:
  * single-mode laws of P Pi_i P for K modes (i = N'-1 is the uplift-decoder mode): free cumulants kappa_2..6,
    rho_n = kappa_n/kappa_n^{Bern(b)}, r = rho_2, exact atom counts;
  * pair invariants (c, F, X) for all pairs among the K modes;
  * 6-letter pair-word ratios for the first `--triples` triples;
  * multi-mode slots Pi_T, |T| = k in --kmodes: single invariants for T = {N'-k..N'-1} and pair invariants
    against the disjoint set {0..k-1};
  * the decoder spectrum (mode N'-1) as .npy.

Output (never overwritten; one JSONL line per realization is appended):
  <out>/raw_q{q}_N{N'}.jsonl, <out>/spectra/..., <out>/meta_q{q}_N{N'}_<timestamp>.json

Example:
  .venv/bin/python scripts/run_chord_invariants.py --q 3 --Np 12 --seeds 0-15 --haar-seeds 0-7 --modes 6 --triples 4 --kmodes 2 3
"""
import argparse
import datetime
import json
import os
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from chord_invariants import (BPSProjector, estimate_memory_gb, pair_invariants, single_invariants,  # noqa: E402
                              triple_invariants)


def parse_range(s):
    out = []
    for part in s.split(","):
        if "-" in part:
            lo, hi = part.split("-")
            out += list(range(int(lo), int(hi) + 1))
        elif part:
            out.append(int(part))
    return out


def git_commit():
    try:
        root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
        h = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
        dirty = subprocess.check_output(["git", "status", "--porcelain"], cwd=root, text=True).strip() != ""
        return h + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


def measure(args, seed, haar):
    t0 = time.time()
    pr = BPSProjector(args.Np, args.P, q=args.q, seed=seed, method=args.method, haar=haar,
                      null_method=args.null_method, cache_gb=args.cache_gb)
    Np = args.Np
    modes = [Np - 1] + list(range(0, args.modes - 1))
    xs = {i: pr.slot([i]) for i in modes}
    row = dict(q=args.q, Np=Np, P=args.P, seed=seed, model="haar" if haar else "syk", method=pr.method,
               null_method=pr.null_method,
               D=pr.D, d=pr.d, r_rank=pr.r, a=pr.a, b=1 - args.P / Np, single=[], pairs=[], triples=[],
               multimode={})
    for i in modes:
        s, lam = single_invariants(pr, xs[i])
        s["mode"] = i
        row["single"].append(s)
        if i == Np - 1 and args.save_spectra:
            sd = os.path.join(args.out, "spectra")
            os.makedirs(sd, exist_ok=True)
            np.save(os.path.join(sd, f"q{args.q}_N{Np}_{row['model']}_s{seed}.npy"), lam)
    for a_i in range(len(modes)):
        for b_i in range(a_i + 1, len(modes)):
            i, j = modes[a_i], modes[b_i]
            p = pair_invariants(pr, xs[i], xs[j])
            p.update(i=i, j=j)
            row["pairs"].append(p)
    trip = [(modes[0], modes[1], modes[2])] + [tuple(modes[t:t + 3]) for t in range(1, len(modes) - 2)]
    for (i, j, k) in trip[: args.triples]:
        row["triples"].append(dict(i=i, j=j, k=k, words=triple_invariants(pr, xs[i], xs[j], xs[k])))
    for k in args.kmodes:
        if 2 * k > Np or args.P > Np - k:
            continue
        T1, T2 = list(range(Np - k, Np)), list(range(0, k))
        y1, y2 = pr.slot(T1), pr.slot(T2)
        s, _ = single_invariants(pr, y1)
        p = pair_invariants(pr, y1, y2)
        row["multimode"][str(k)] = dict(T1=T1, T2=T2, single=s, pair=p)
    row["secs"] = time.time() - t0
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--q", type=int, default=3)
    ap.add_argument("--Np", type=int, required=True, help="enlarged size N'")
    ap.add_argument("--P", type=int, default=None, help="charge (default floor(N'/2))")
    ap.add_argument("--seeds", default="0-3", help="SYK coupling seeds, e.g. 0-15 or 0,2,5")
    ap.add_argument("--haar-seeds", default="", help="seeds for Haar-projector null model")
    ap.add_argument("--modes", type=int, default=6, help="number of single modes (includes N'-1)")
    ap.add_argument("--triples", type=int, default=4)
    ap.add_argument("--kmodes", type=int, nargs="*", default=[2, 3])
    ap.add_argument("--method", default="auto", choices=["auto", "basis", "complement"])
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                                                  "results", "data", "chord_invariants_2026-09-28"))
    ap.add_argument("--no-spectra", dest="save_spectra", action="store_false")
    ap.add_argument("--null-method", default="auto", choices=["auto", "svd", "eigh"],
                    help="BPS basis from full SVD of [Q;Q^dag] or eigh of {Q,Q^dag} (auto: svd iff D <= 7000)")
    ap.add_argument("--cache-gb", type=float, default=4.0, help="byte budget of the compressed-letter LRU cache")
    ap.add_argument("--dry-run", action="store_true", help="print the peak-memory estimate and exit")
    args = ap.parse_args()
    args.P = args.Np // 2 if args.P is None else args.P
    meth = args.method if args.method != "auto" else ("complement" if args.q >= 5 else "basis")
    est = estimate_memory_gb(args.Np, args.P, args.q, method=meth, null_method=args.null_method,
                             cache_gb=args.cache_gb)
    print(f"memory estimate (run under scripts/memwatch.py with a cap above this): {est}", flush=True)
    if args.dry_run:
        return
    args.out = os.path.abspath(args.out)
    os.makedirs(args.out, exist_ok=True)
    stamp = datetime.datetime.now().strftime("%Y%m%dT%H%M%S")
    meta = dict(vars(args), git=git_commit(), started=stamp, host_numpy=np.__version__,
                description="chord invariants of BPS-projected slot projectors; see src/chord_invariants.py")
    with open(os.path.join(args.out, f"meta_q{args.q}_N{args.Np}_{stamp}.json"), "w") as f:
        json.dump(meta, f, indent=1)
    raw = os.path.join(args.out, f"raw_q{args.q}_N{args.Np}.jsonl")
    jobs = [(s, False) for s in parse_range(args.seeds)] + [(s, True) for s in parse_range(args.haar_seeds)]
    for seed, haar in jobs:
        row = measure(args, seed, haar)
        with open(raw, "a") as f:
            f.write(json.dumps(row) + "\n")
        s0 = row["single"][0]
        X = np.mean([p["X"] for p in row["pairs"]]) if row["pairs"] else float("nan")
        print(f"q={args.q} N'={args.Np} {row['model']} seed={seed}: a={row['a']:.3f} r={s0['r']:.4f} "
              f"rho4={s0['rho'][3]:.4f} X={X:.4f} ({row['secs']:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
