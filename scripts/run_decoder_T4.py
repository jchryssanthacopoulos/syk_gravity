#!/usr/bin/env python3
r"""
Decoder fourth moment T_4 = tau((P U P)^4) vs its exact size decomposition (research note section 11):
    T_4(U) = c_0 + 2 a (T_2(U) - a) + R(U)                      (exact, per realization, any U in U(N))
    c_0 = sum_k phi_k w_k  (two-independent-model value),  R(U) = <U P_{>=1} U, delta_P U P_{>=1} U>/d
and the size-diagonal approximation  R(U) ~ sum_{k,k'>=1} chi_hat_k(U) chi_hat_k'(U) Gamma_kk'.
Probes: single-site parities U_i for several sites i (the decoder) and parity sets U_S, S = {N'-1, 0, .., s-2}.
Output: results/data/parity_family_2026-09-29/decoderT4_q{q}_N{N'}.jsonl (append).
Memory ~ (P+4) dense D x D complex matrices (N'=14: ~2.5 GB).  Run under scripts/memwatch.py.
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
from size_decomposition import (HopTable, chi_hat, decoder_T4_pieces, size_components,  # noqa: E402
                                transmission_matrix, transmissions_from_weights)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--q", type=int, default=3)
    ap.add_argument("--Np", type=int, required=True)
    ap.add_argument("--seeds", default="0-3")
    ap.add_argument("--haar-seeds", default="")
    ap.add_argument("--sites", type=int, default=4, help="number of single sites for the decoder (s = 1)")
    ap.add_argument("--out", default=os.path.join(ROOT, "results/data/parity_family_2026-09-29"))
    a = ap.parse_args()

    def rng(s):
        if not s:
            return []
        lo, _, hi = s.partition("-")
        return list(range(int(lo), int(hi or lo) + 1))
    N, P = a.Np, a.Np // 2
    table = HopTable(N, P)
    A = transmission_matrix(N, P, table.masks)
    for seed, haar in [(s, False) for s in rng(a.seeds)] + [(s, True) for s in rng(a.haar_seeds)]:
        t0 = time.time()
        pr = BPSProjector(N, P, q=a.q, seed=seed, method="basis", haar=haar)
        Pm = pr.B @ pr.B.conj().T
        comps, resid = size_components(Pm, table)
        d = pr.d
        w0 = np.array([float(np.real(np.vdot(C, C))) / d for C in comps])
        phi = transmissions_from_weights(w0, d, N, P, A)
        c0, G, w = decoder_T4_pieces(Pm, comps, phi)
        del comps
        R1 = float(G.sum())
        probes = [("site", [i]) for i in ([N - 1] + list(range(a.sites - 1)))]
        probes += [("parity", [N - 1] + list(range(s - 1))) for s in range(2, N // 2 + 1)]
        res = []
        for kind, S in probes:
            s = len(S)
            u = 2 * parity_projector_diag(pr.masks, S) - 1
            M = pr.B.conj().T @ (u[:, None] * pr.B)
            mu = np.linalg.eigvalsh((M + M.conj().T) / 2)
            T2, T4 = float(np.mean(mu ** 2)), float(np.mean(mu ** 4))
            ch = np.array([chi_hat(N, s, k) for k in range(len(w))])
            T2pred = float(np.dot(w, ch))
            Rexact = T4 - c0 - 2 * pr.a * (T2 - pr.a)
            Rdiag = float(ch[1:] @ G @ ch[1:])
            res.append(dict(kind=kind, S=S, s=s, T2=T2, T4=T4, T2pred=T2pred, R_exact=Rexact, R_diag=Rdiag,
                            T4_pred_diag=c0 + 2 * pr.a * (T2pred - pr.a) + Rdiag, chi=ch.tolist()))
        row = dict(q=a.q, Np=N, P=P, seed=seed, model="haar" if haar else "syk", a=pr.a, d=d, w=w.tolist(),
                   phi=phi.tolist(), c0=c0, Gamma=G.tolist(), R_identity=R1,
                   R_identity_check=1 - c0 - 2 * pr.a * (1 - pr.a), resid=resid, probes=res, secs=time.time() - t0)
        with open(os.path.join(a.out, f"decoderT4_q{a.q}_N{N}.jsonl"), "a") as f:
            f.write(json.dumps(row) + "\n")
        print(f"N'={N} {row['model']} seed={seed}: c0={c0:.5f}  R(1)={R1:.5f} vs 1-c0-2a(1-a)={row['R_identity_check']:.5f}"
              f"  ({row['secs']:.0f}s)", flush=True)
        for e in res:
            print(f"   {e['kind']:6s} s={e['s']} S={e['S'][:3]}..: T4={e['T4']:.5f}  R_exact={e['R_exact']:+.5f}  R_diag={e['R_diag']:+.5f}"
                  f"  T4_pred={e['T4_pred_diag']:.5f}", flush=True)


if __name__ == "__main__":
    main()
