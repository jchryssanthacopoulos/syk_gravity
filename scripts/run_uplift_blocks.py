#!/usr/bin/env python3
r"""
First order -> self-dressed -> exact decomposition of the decoder via the block form (docs/derivations.md D6 §5).

On B = B_h (+) B_g (exact orthogonal splitting) the decoder mu-operator is [[M_h, C2], [C2^dag, M_g]].
  first order:   spectrum (1 - sigma^2)/(1 + sigma^2) and -(1 - sigma'^2)/(1 + sigma'^2)   (C2 = 0, no self-dressing)
  self-dressed:  spec M_h  U  spec M_g   (exact uplift within each family, h-g coupling dropped)
  exact:         spec of the full block operator (identical to the N'-mode decoder; tested in tests/test_uplift.py)
Reports T_2 = <y>, H = <y(1-y)>, T_4 = <y^2> (y = mu^2), frac(y < 1/2), and the mixing term 2||C2||_F^2 / d.

Usage: .venv/bin/python scripts/run_uplift_blocks.py --Np 10 11 12 13 14 15 --nreal 2 --seed0 9000
Memory: old-model dense maps only (N' = 15: ~1 GB; N' = 16: ~3 GB).
"""
import argparse
import json
import os
import sys

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from chord_chaos import Sector, random_couplings  # noqa: E402
from uplift import OldModel  # noqa: E402
from run_uplift_walls import git_state  # noqa: E402


def stats(mu):
    y = np.asarray(mu) ** 2
    return dict(T2=float(y.mean()), H=float(np.mean(y * (1 - y))), T4=float(np.mean(y * y)), frac_low=float(np.mean(y < 0.5)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Np", type=int, nargs="+", required=True)
    ap.add_argument("--nreal", type=int, default=2)
    ap.add_argument("--seed0", type=int, default=9000)
    ap.add_argument("--out", default=os.path.join(ROOT, "results/data/uplift_2026-10-12/blocks.jsonl"))
    a = ap.parse_args()
    for Np in a.Np:
        for r in range(a.nreal):
            seed = a.seed0 + r
            rng = np.random.default_rng(seed)
            c = random_couplings(rng, Sector(Np, 3, F=Np // 2).M)
            om = OldModel(Np, c, rng)
            sig, sigp = om.first_order_singular_values()
            mu1 = np.concatenate([(1 - sig ** 2) / (1 + sig ** 2), -(1 - sigp ** 2) / (1 + sigp ** 2)])
            Mh, Mg, C2, overlap = om.decoder_blocks()
            mh, mg = np.linalg.eigvalsh(Mh), np.linalg.eigvalsh(Mg)
            mu = np.linalg.eigvalsh(np.block([[Mh, C2], [C2.conj().T, Mg]]))
            d = len(mu)
            rec = dict(Np=Np, seed=seed, d=d, d_h=len(mh), d_g=len(mg), overlap=overlap,
                       first=stats(mu1), self=stats(np.concatenate([mh, mg])), exact=stats(mu),
                       mixing_T2=float(2 * np.linalg.norm(C2) ** 2 / d),
                       mu_h_range=[float(mh.min()), float(mh.max())], mu_g_range=[float(mg.min()), float(mg.max())],
                       git=git_state())
            with open(a.out, "a") as fh:
                fh.write(json.dumps(rec) + "\n")
            f, s, e = rec["first"], rec["self"], rec["exact"]
            print(f"N'={Np} seed {seed}: T2 {f['T2']:.4f} -> {s['T2']:.4f} -> {e['T2']:.4f} (mixing {rec['mixing_T2']:.4f});"
                  f"  H {f['H']:.4f} -> {s['H']:.4f} -> {e['H']:.4f};  |<h|g>| {overlap:.1e}", flush=True)


if __name__ == "__main__":
    main()
