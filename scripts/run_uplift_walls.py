#!/usr/bin/env python3
r"""
Decoder walls from the uplift map (docs/derivations.md D6).

For each realization of one-flavor N=2 SYK at N' modes (p = 3, F = N'//2, added mode = last mode):
  - first-order uplift prediction from the OLD model (N'-1 modes) alone: lambda = 1/(1+sigma^2) for sigma in
    sv(K|B_o^F), K = Q0^+ eta^ ; lambda = sigma'^2/(1+sigma'^2) for sigma' in sv(K'|B_o^{F-1});
  - (optional --exact) the exact decoder spectrum of the N'-mode model, and per eigenvector the uplift obstruction
    r(v) = ||eta ^ P_B a|| / ||P_B a||  (a = n_e = 0 component, P_B = projector on B_o^F);
  - wall fractions (min(lambda, 1-lambda) < eps), exact atom counts, and the refined atom bound / dimension ratio.

Usage: .venv/bin/python scripts/run_uplift_walls.py --Np 14 --nreal 3 --exact --out results/data/uplift_2026-10-12
Memory: old-model dense maps C(N'-1, F+2) x C(N'-1, F-1) etc.; exact needs the N'-mode BPS basis (N'=15: ~2 GB,
N'=16: ~7 GB).  Run heavy cases under scripts/memwatch.py.
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
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from chord_chaos import Sector, random_couplings  # noqa: E402
from run_ou_fourarc import bps_basis_fast  # noqa: E402
from uplift import OldModel, atom_bound, bps_dim, decoder_from_basis, wall_fractions  # noqa: E402

EPS = (0.1, 0.03, 0.01, 0.001)


def git_state():
    try:
        g = subprocess.run(["git", "-C", ROOT, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", ROOT, "status", "--porcelain"], capture_output=True, text=True).stdout.strip()
        return g + ("-dirty" if dirty else "")
    except OSError:
        return "unknown"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Np", type=int, required=True)
    ap.add_argument("--nreal", type=int, default=3)
    ap.add_argument("--seed0", type=int, default=9000, help="realization r uses numpy seed seed0 + r")
    ap.add_argument("--exact", action="store_true", help="also diagonalize the N'-mode model")
    ap.add_argument("--out", default=os.path.join(ROOT, "results/data/uplift_2026-10-12"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    p, F = 3, a.Np // 2
    SN = Sector(a.Np, p, F=F)
    bound = atom_bound(a.Np, p, F)
    path = os.path.join(a.out, f"walls_N{a.Np}.jsonl")
    meta = dict(args=vars(a), git=git_state(), atom_bound=bound, eps=EPS)
    with open(os.path.join(a.out, f"meta_walls_N{a.Np}_{time.strftime('%Y%m%dT%H%M%S')}.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    for r in range(a.nreal):
        t0 = time.time()
        rng = np.random.default_rng(a.seed0 + r)
        c = random_couplings(rng, SN.M)
        om = OldModel(a.Np, c, rng)
        lam1 = om.first_order_decoder()
        rec = dict(Np=a.Np, real=r, seed=a.seed0 + r, first=wall_fractions(lam1, EPS), first_lam=np.sort(lam1).tolist(),
                   first_atoms=int(np.sum(np.minimum(lam1, 1 - lam1) < 1e-9)))
        if a.exact:
            V, _, _ = bps_basis_fast(SN, c, bps_dim(a.Np, p, F), rng)
            u = 1.0 - 2.0 * np.array([(m >> (a.Np - 1)) & 1 for m in SN.masks], float)
            mu, W = np.linalg.eigh(V.conj().T @ (u[:, None] * V))
            lam = (1 + mu) / 2
            X = V @ W
            # uplift obstruction per eigenvector: harmonic part of the n_e = 0 component
            idx = {m: k for k, m in enumerate(om.SF.masks)}
            rows0 = [k for k, m in enumerate(SN.masks) if not (m >> (a.Np - 1)) & 1]
            A = np.zeros((om.SF.D, X.shape[1]), complex)
            A[[idx[SN.masks[k]] for k in rows0], :] = X[rows0, :]
            H = om.VF @ (om.VF.conj().T @ A)
            nh = np.linalg.norm(H, axis=0)
            ratio = np.linalg.norm(om.E_F @ H, axis=0) / np.maximum(nh, 1e-300)
            rec.update(exact=wall_fractions(lam, EPS), lam=lam.tolist(), h_weight=(nh ** 2).tolist(),
                       obstruction=ratio.tolist(), exact_atoms=int(np.sum(np.minimum(lam, 1 - lam) < 1e-9)),
                       T=[float(np.mean((2 * lam - 1) ** (2 * k))) for k in range(1, 5)])
        rec["secs"] = time.time() - t0
        with open(path, "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        print(a.Np, r, f"{rec['secs']:.0f}s", "first", np.round(rec["first"], 4), "exact",
              np.round(rec.get("exact", [np.nan]), 4), "atoms", rec.get("exact_atoms"), flush=True)


if __name__ == "__main__":
    main()
