#!/usr/bin/env python3
r"""
Out-of-sample test of the LMRS R-charge dependence (docs/derivations.md D1).

For even N' the BPS states sit in three fermion-number sectors: P = N'/2 (j = 0) and P = N'/2 +- 1 (j = +-1/3).
The zero-energy variance of the single-mode decoder in sector j is predicted (no free parameter) to be
    Var_j / Var_0 = cos(pi/3) Gamma(7/6) Gamma(1/2) / Gamma(5/6)^2 = 0.645      (LMRS formula)
                                                                    ~ 0.66      (UV-capped toy, var_capped)
PRE-REGISTERED 2026-09-28, before any j = +-1/3 measurement (see the D1 note for the table by N').
Same couplings (seeds) are used in all three sectors, so the ratio is paired.

Usage: .venv/bin/python scripts/memwatch.py --limit-gb 3 -- .venv/bin/python scripts/lmrs_jsector_test.py --Np 8,10,12,14
Output: results/data/lmrs_variance_2026-09-28/jsector.jsonl (append).
"""
import argparse
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from chord_invariants import BPSProjector  # noqa: E402
from lmrs_predictions import var_capped, var_zero_energy  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Np", default="8,10,12,14")
    ap.add_argument("--seeds", default="0-3")
    ap.add_argument("--modes", type=int, default=4)
    ap.add_argument("--out", default=os.path.join(ROOT, "results/data/lmrs_variance_2026-09-28/jsector.jsonl"))
    a = ap.parse_args()
    lo, _, hi = a.seeds.partition("-")
    seeds = range(int(lo), int(hi or lo) + 1)
    for Np in (int(x) for x in a.Np.split(",")):
        assert Np % 2 == 0, "j = +-1/3 sectors are for even N'"
        for s in seeds:
            t0 = time.time()
            row = dict(Np=Np, seed=s, sectors={})
            for P in (Np // 2, Np // 2 - 1, Np // 2 + 1):
                pr = BPSProjector(Np, P, q=3, seed=s, method="basis")
                v = []
                for i in [Np - 1] + list(range(a.modes - 1)):
                    x = pr.slot([i])
                    A = pr.B.conj().T @ (x[:, None] * pr.B)
                    lam = np.linalg.eigvalsh((A + A.conj().T) / 2)
                    v.append(float(lam.var()))
                b = 1 - P / Np
                row["sectors"][str(P)] = dict(j=(P - Np / 2) / 3, d=pr.d, D=pr.D, var=float(np.mean(v)),
                                              var_modes=v, b=b, pred_lmrs=var_zero_energy(Np, P),
                                              pred_capped=var_capped(Np, b * (1 - b), P))
            S = row["sectors"]
            row["ratio_minus"] = S[str(Np // 2 - 1)]["var"] / S[str(Np // 2)]["var"]
            row["ratio_plus"] = S[str(Np // 2 + 1)]["var"] / S[str(Np // 2)]["var"]
            row["secs"] = time.time() - t0
            with open(a.out, "a") as f:
                f.write(json.dumps(row) + "\n")
            print(f"N'={Np} seed={s}: d(j=0,-1/3,+1/3)=" + ",".join(str(S[k]["d"]) for k in S) +
                  f"  Var_j/Var_0: -1/3 {row['ratio_minus']:.3f}, +1/3 {row['ratio_plus']:.3f}"
                  f"  (pred 0.645 / toy ~0.66)  [{row['secs']:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
