#!/usr/bin/env python3
r"""
Size of the crossed matter-chord diagram in the third moment of the first-order dressing operator (docs/derivations.md D7).

E_eta tr D^3 is a sum over the 3! Wick pairings of (E, E^dag); five are non-crossing, one is crossing.  With three
independent Gaussian vectors a, b, g (E|.|^2 = 2), X_v = G^1/2 E_v V_b, the estimator
    (1/b) Tr[(X_g^dag X_a)(X_b^dag X_g)(X_a^dag X_b)]
is unbiased for the crossing term alone (each vector is contracted only with itself).  The identity and full-cycle
pairings are estimated the same way and compared with their exact values (tr Dbar^3 and (1/b) Tr S^3) as a check.

Usage: .venv/bin/python scripts/run_chord_crossing.py --Np 12 --nsamp 1000 --seeds 9000 9001
"""
import argparse
import json
import os
import sys
import time
from math import comb

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from chord_dressing import DressingFamily, crossing_terms, old_model_lean  # noqa: E402
from run_chord_dressing import OUT, git_state  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Np", type=int, required=True)
    ap.add_argument("--seeds", type=int, nargs="+", default=[9000])
    ap.add_argument("--nsamp", type=int, default=500)
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, f"crossing_N{a.Np}.jsonl")
    with open(os.path.join(a.out, f"meta_crossing_N{a.Np}_{time.strftime('%Y%m%dT%H%M%S')}.json"), "w") as fh:
        json.dump(dict(args=vars(a), git=git_state()), fh, indent=1)
    for seed in a.seeds:
        om = old_model_lean(a.Np, seed)
        rng = np.random.default_rng(seed + 99)
        rec = dict(Np=a.Np, seed=seed, nsamp=a.nsamp)
        for fam in ("h", "g"):
            Fm = DressingFamily(om, fam)
            V = Fm.V_inner()
            Vt = Fm.Ue.conj().T @ V @ Fm.Ue
            sev = np.linalg.eigvalsh(((Vt + Vt.conj().T) / 2) / np.sqrt(np.outer(Fm.we, Fm.we)))
            Dbar = Fm.lift(Fm.G_inner())
            r = crossing_terms(Fm, a.nsamp, rng, comb(om["N"], 2))
            r["identity_exact"] = float(np.trace(Dbar @ Dbar @ Dbar).real / Fm.b)
            r["cycle_exact"] = float(np.sum(sev ** 3) / Fm.b)
            rec[fam] = r
            print(f"N'={a.Np} seed={seed} {fam}: crossing {r['crossing'][0]:.3e} +- {r['crossing'][1]:.1e} "
                  f"({100 * r['crossing'][0] / r['m3'][0]:+.2f}% of m3 = {r['m3'][0]:.4g}) | identity {r['identity'][0]:.4g} "
                  f"(exact {r['identity_exact']:.4g}) | cycle {r['cycle'][0]:.4g} (exact {r['cycle_exact']:.4g})", flush=True)
        with open(path, "a") as fh:
            fh.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main()
