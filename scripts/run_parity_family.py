#!/usr/bin/env python3
r"""
Parity-family decoder spectra (research/notes/decoder_moments_2026-09-29.md).

For each realization (SYK couplings seed, or a Haar projector of the same rank) and each s = 1..smax, compute the
full spectrum of P_B Pi_S P_B on the charge-P BPS space, with Pi_S = (1 + (-1)^{N_S})/2 and
S = {N'-1, 0, 1, ..., s-2} (s = 1 is the ordinary uplift decoder).  The per-Q-chord crossing weight of U_S is
x_s (exact: hypergeometric; double-scaled: q^{s/p}); x = 0 exactly at s = N'/2 for even N', odd p.

Records per (realization, s): a, b_S, lambda moments (1..8), mu = 2 lambda - 1 raw moments (1..8), free cumulants
and rho_n (n = 2..6), KS to Wachter(a, b_S), atom counts; spectra saved as .npy.
Output: <out>/raw_q{q}_N{N'}.jsonl (append), <out>/spectra/, <out>/meta_*.json.

Example: .venv/bin/python scripts/memwatch.py --limit-gb 4 -- \
         .venv/bin/python scripts/run_parity_family.py --Np 12 --seeds 0-3 --haar-seeds 0
"""
import argparse
import datetime
import json
import os
import subprocess
import sys
import time

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from chord_invariants import BPSProjector, estimate_memory_gb, single_invariants  # noqa: E402
from nullmodels import ks_wachter  # noqa: E402
from parity_decoder import mu_moments_from_lambda, parity_projector_diag, x_double_scaled, x_exact  # noqa: E402


def parse_range(s):
    out = []
    for part in s.split(","):
        if part:
            lo, _, hi = part.partition("-")
            out += list(range(int(lo), int(hi or lo) + 1))
    return out


def git_commit():
    try:
        h = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
        dirty = subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).strip() != ""
        return h + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--q", type=int, default=3)
    ap.add_argument("--Np", type=int, required=True)
    ap.add_argument("--P", type=int, default=None)
    ap.add_argument("--seeds", default="0-3")
    ap.add_argument("--haar-seeds", default="")
    ap.add_argument("--smax", type=int, default=None, help="default N'//2")
    ap.add_argument("--method", default="basis", choices=["basis", "complement", "auto"])
    ap.add_argument("--out", default=os.path.join(ROOT, "results/data/parity_family_2026-09-29"))
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    a.P = a.Np // 2 if a.P is None else a.P
    a.smax = a.Np // 2 if a.smax is None else a.smax
    est = estimate_memory_gb(a.Np, a.P, a.q, method=a.method if a.method != "auto" else "basis")
    print(f"memory estimate: {est}", flush=True)
    if a.dry_run:
        return
    os.makedirs(os.path.join(a.out, "spectra"), exist_ok=True)
    stamp = datetime.datetime.now().strftime("%Y%m%dT%H%M%S")
    with open(os.path.join(a.out, f"meta_q{a.q}_N{a.Np}_{stamp}.json"), "w") as f:
        json.dump(dict(vars(a), git=git_commit(), started=stamp,
                       description="parity-family decoder spectra; see scripts/run_parity_family.py"), f, indent=1)
    jobs = [(s, False) for s in parse_range(a.seeds)] + [(s, True) for s in parse_range(a.haar_seeds)]
    for seed, haar in jobs:
        t0 = time.time()
        pr = BPSProjector(a.Np, a.P, q=a.q, seed=seed, method=a.method, haar=haar)
        model = "haar" if haar else "syk"
        row = dict(q=a.q, Np=a.Np, P=a.P, seed=seed, model=model, D=pr.D, d=pr.d, a=pr.a, family=[])
        for s in range(1, a.smax + 1):
            S = [a.Np - 1] + list(range(s - 1))
            xd = parity_projector_diag(pr.masks, S)
            inv, lam = single_invariants(pr, xd, K=6)
            b = inv["b_uv"]
            mom = [float(np.mean(lam ** k)) for k in range(1, 9)]
            row["family"].append(dict(
                s=s, S=S, b=b, x_exact=x_exact(a.Np, a.q, s), x_ds=x_double_scaled(a.Np, a.q, s),
                lam_moments=mom, mu_moments=mu_moments_from_lambda(lam, 8), kappa=inv["kappa"], rho=inv["rho"],
                r=inv["r"], n0=inv["n0"], n1=inv["n1"], ks_wachter=ks_wachter(lam, pr.a, b)))
            np.save(os.path.join(a.out, "spectra", f"q{a.q}_N{a.Np}_{model}_seed{seed}_s{s}.npy"), lam)
        row["secs"] = time.time() - t0
        with open(os.path.join(a.out, f"raw_q{a.q}_N{a.Np}.jsonl"), "a") as f:
            f.write(json.dumps(row) + "\n")
        fam = row["family"]
        print(f"q={a.q} N'={a.Np} {model} seed={seed} a={pr.a:.3f} ({row['secs']:.0f}s)  s: r / KS_W = " +
              "  ".join(f"{e['s']}:{e['r']:.3f}/{e['ks_wachter']:.3f}" for e in fam), flush=True)


if __name__ == "__main__":
    main()
