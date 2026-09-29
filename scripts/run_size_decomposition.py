#!/usr/bin/env python3
r"""
Operator-size (U(N)-irrep) decomposition of the BPS projector and the exact finite-N sum rule for T_2
(research/notes/decoder_moments_2026-09-29.md).

For each realization: P = B B^dag on Lambda^P; w_k = ||P^{(k)}||^2 / d (k = 0..P; w_0 = a exactly).
Prediction (U(N)-averaged, exact in expectation):  E[T_2(s)] = sum_k w_k chi_hat_k(U_S),  T_2 = tau(P U_S P U_S).
Also records the measured T_2(s) = tau(P U_S P U_S) for the same realization and S = {N'-1, 0, .., s-2}.

Output: results/data/parity_family_2026-09-29/size_q{q}_N{N'}.jsonl (append).
Memory: dense D x D complex matrices, about (P+4) copies: N'=14 -> ~2 GB.  Run under scripts/memwatch.py.
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
from parity_decoder import parity_projector_diag  # noqa: E402
from size_decomposition import HopTable, chi_hat, size_weights  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--q", type=int, default=3)
    ap.add_argument("--Np", type=int, required=True)
    ap.add_argument("--seeds", default="0-3")
    ap.add_argument("--haar-seeds", default="")
    ap.add_argument("--out", default=os.path.join(ROOT, "results/data/parity_family_2026-09-29"))
    a = ap.parse_args()

    def rng(s):
        if not s:
            return []
        lo, _, hi = s.partition("-")
        return list(range(int(lo), int(hi or lo) + 1))
    P = a.Np // 2
    table = HopTable(a.Np, P)
    for seed, haar in [(s, False) for s in rng(a.seeds)] + [(s, True) for s in rng(a.haar_seeds)]:
        t0 = time.time()
        pr = BPSProjector(a.Np, P, q=a.q, seed=seed, method="basis", haar=haar)
        Pm = pr.B @ pr.B.conj().T
        w, resid = size_weights(Pm, table)
        fam = []
        for s in range(1, a.Np // 2 + 1):
            S = [a.Np - 1] + list(range(s - 1))
            u = 2 * parity_projector_diag(pr.masks, S) - 1          # U_S diagonal (+-1)
            A = pr.B.conj().T @ (u[:, None] * pr.B)
            T2 = float(np.real(np.trace(A @ A))) / pr.d
            pred = float(sum(w[k] * chi_hat(a.Np, s, k) for k in range(len(w))))
            fam.append(dict(s=s, T2_measured=T2, T2_pred=pred, chi=[chi_hat(a.Np, s, k) for k in range(len(w))]))
        row = dict(q=a.q, Np=a.Np, P=P, seed=seed, model="haar" if haar else "syk", a=pr.a, d=pr.d,
                   w=w.tolist(), resid=resid, family=fam, secs=time.time() - t0)
        with open(os.path.join(a.out, f"size_q{a.q}_N{a.Np}.jsonl"), "a") as f:
            f.write(json.dumps(row) + "\n")
        print(f"N'={a.Np} {row['model']} seed={seed} a={pr.a:.4f} w0={w[0]:.4f} w=" + " ".join(f"{x:.4f}" for x in w) +
              f" (resid {resid:.1e}; {row['secs']:.0f}s)", flush=True)
        print("   s: T2 measured / predicted  " + "  ".join(f"{e['s']}:{e['T2_measured']:.4f}/{e['T2_pred']:.4f}" for e in fam),
              flush=True)


if __name__ == "__main__":
    main()
