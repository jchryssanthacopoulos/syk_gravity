#!/usr/bin/env python3
r"""
Diagnose why exact SYK decoder variances undershoot the parameter-free LMRS prediction (docs/derivations.md D1).

Per realization (q_hat = 3, charge P = floor(N'/2), same couplings/seeds as the campaign):
  1. BPS gap -> effective Schwarzian coupling.  LMRS eq. (81): 8 C E_gap = (j_c - 1/2)^2, minimized over the
     continuum multiplets reaching our R-charge J (J = 0: 1/4;  J = -1/6: 1/9, via L_{s,j'} with j'-1 = J).
     Couplings: our C_T are unit complex Gaussians, E|C|^2 = 2 = 2J/N'^2 (FGMS eq. 5.4; LMRS eq. 77 prints 2J/N,
     a typo: it would make H non-extensive)  =>  J = N'^2.  So alpha_S,eff = C_eff J / N' = C_eff N'.
  2. Zero-energy two-point functions in the BPS sector of
       same-site  O = n_i - mu             (our decoder variance)
       hopping    A = c_k^dag c_i, i != k  (the operator LMRS eq. 85 actually treats)
     LMRS predict the same value for both at leading large N.
Outputs one JSON line per realization to --out (append).  Dense linear algebra, D <= C(15,7) = 6435 (< 1 GB).

Usage: .venv/bin/python scripts/memwatch.py --limit-gb 4 -- .venv/bin/python scripts/lmrs_mismatch_checks.py \
          --Np 8-15 --seeds 0-3
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from q_scan import form_wedge, rand_form, subset_index  # noqa: E402


def parse_range(s):
    out = []
    for part in s.split(","):
        lo, _, hi = part.partition("-")
        out += list(range(int(lo), int(hi or lo) + 1))
    return out


def hop(masks, index, i, k):
    """Sparse matrix of c_k^dag c_i on the fixed-charge basis (bitmask ordering as subset_index)."""
    I, J, V = [], [], []
    for col, S in enumerate(masks):
        if not (S >> i) & 1 or (S >> k) & 1:
            continue
        s1 = bin(S & ((1 << i) - 1)).count("1")
        S2 = S & ~(1 << i)
        s2 = bin(S2 & ((1 << k) - 1)).count("1")
        I.append(index[S2 | (1 << k)])
        J.append(col)
        V.append((-1.0) ** (s1 + s2))
    return sp.csr_matrix((V, (I, J)), shape=(len(masks), len(masks)))


def run(Np, seed, q=3):
    t0 = time.time()
    P = Np // 2
    masks, index = subset_index(Np, P)
    D = len(masks)
    C = rand_form(seed, Np, q)
    blocks = [form_wedge(C, Np, P, q)]
    if P - q >= 0:
        blocks.append(form_wedge(C, Np, P - q, q).conj().T)
    Ms = sp.vstack(blocks).tocsr()
    H = (Ms.conj().T @ Ms).toarray()
    w, V = np.linalg.eigh(H)
    tol = D * w[-1] * np.finfo(float).eps * 100
    d = int((w <= tol).sum())
    E_gap = float(w[d])
    B = V[:, :d]
    j = (P - Np / 2) / q
    num = 0.25 if j == 0 else min((j - 0.5) ** 2, (j + 1 - 0.5) ** 2)     # multiplets H_{j}, L_{j+1}
    C_eff = num / (8 * E_gap)
    J = Np ** 2                                                            # E|C|^2 = 2 = 2J/N'^2
    out = dict(q=q, Np=Np, P=P, seed=seed, D=D, d=d, j=j, E_gap=E_gap, C_eff=C_eff,
               alpha_S_eff=C_eff * J / Np, same=[], hop=[])
    occ = lambda i: np.array([(m >> i) & 1 for m in masks], float)         # noqa: E731
    for i in (Np - 1, 0, 1):
        x = occ(i)
        A = B.conj().T @ (x[:, None] * B)
        mu = np.trace(A).real / d
        out["same"].append(float(np.trace(A @ A).real / d - mu ** 2))
    for (i, k) in ((Np - 1, 0), (0, 1), (1, 2)):
        Hk = hop(masks, index, i, k)
        A = B.conj().T @ (Hk @ B)
        m1 = np.trace(A) / d
        out["hop"].append(float(np.trace(A @ A.conj().T).real / d - abs(m1) ** 2))
    out["secs"] = time.time() - t0
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Np", default="8-14")
    ap.add_argument("--seeds", default="0-3")
    ap.add_argument("--out", default=os.path.join(ROOT, "results/data/lmrs_variance_2026-09-28/mismatch_checks.jsonl"))
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    for Np in parse_range(a.Np):
        for s in parse_range(a.seeds):
            r = run(Np, s)
            with open(a.out, "a") as f:
                f.write(json.dumps(r) + "\n")
            print(f"N'={Np} seed={s} d={r['d']} E_gap={r['E_gap']:.4f} alpha_S,eff={r['alpha_S_eff']:.5f} "
                  f"Var_same={np.mean(r['same']):.4f} Var_hop={np.mean(r['hop']):.4f} ({r['secs']:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
