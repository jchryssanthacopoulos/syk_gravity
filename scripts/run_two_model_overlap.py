#!/usr/bin/env python3
r"""
Principal angles between the BPS spaces of two independent one-flavor N=2 SYK models (same N', same charge P).

Motivation (research/notes/decoder_moments_2026-09-29.md): (P U P)^2 = P P' P with P' = U P U the BPS projector of
Q' = U Q U, whose couplings are correlated with those of Q with coefficient x per chord.  At x = 0 (parity set
|S| = N'/2) the double-scaled chord theory treats Q and Q' as independent.  Prediction: the moments
E[nu^k] = tau((P1 P2 P1)^k) of two independent models equal T_{2k} of the x = 0 parity decoder.
Free compression (Haar) gives nu ~ Wachter(a, a).

Output: <out>/two_model_q{q}_N{N'}.jsonl (append): per seed pair, nu moments k = 1..4, KS to Wachter(a, a),
spectrum saved in <out>/spectra/.
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
from nullmodels import ks_wachter  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--q", type=int, default=3)
    ap.add_argument("--Np", type=int, required=True)
    ap.add_argument("--pairs", default="100:101,102:103,104:105,106:107",
                    help="seed pairs (disjoint from the parity-family seeds 0-3)")
    ap.add_argument("--haar", action="store_true", help="also do one Haar-Haar pair (control)")
    ap.add_argument("--out", default=os.path.join(ROOT, "results/data/parity_family_2026-09-29"))
    a = ap.parse_args()
    os.makedirs(os.path.join(a.out, "spectra"), exist_ok=True)
    P = a.Np // 2
    jobs = [(tuple(int(x) for x in pr.split(":")), False) for pr in a.pairs.split(",")]
    if a.haar:
        jobs.append(((200, 201), True))
    for (s1, s2), haar in jobs:
        t0 = time.time()
        p1 = BPSProjector(a.Np, P, q=a.q, seed=s1, method="basis", haar=haar, haar_seed=s1 if haar else None)
        p2 = BPSProjector(a.Np, P, q=a.q, seed=s2, method="basis", haar=haar, haar_seed=s2 if haar else None)
        Mx = p1.B.conj().T @ p2.B                                   # d1 x d2
        nu = np.clip(np.linalg.eigvalsh(Mx @ Mx.conj().T).real, 0, 1)   # cos^2 principal angles (d1 values)
        aa = p1.a
        row = dict(q=a.q, Np=a.Np, P=P, seeds=[s1, s2], model="haar" if haar else "syk", a=aa, d=p1.d,
                   nu_moments=[float(np.mean(nu ** k)) for k in range(1, 5)], ks_wachter_aa=ks_wachter(nu, aa, aa),
                   n0=int((nu < 1e-9).sum()), n1=int((nu > 1 - 1e-9).sum()), secs=time.time() - t0)
        np.save(os.path.join(a.out, "spectra", f"twomodel_q{a.q}_N{a.Np}_{row['model']}_{s1}_{s2}.npy"), nu)
        with open(os.path.join(a.out, f"two_model_q{a.q}_N{a.Np}.jsonl"), "a") as f:
            f.write(json.dumps(row) + "\n")
        m = row["nu_moments"]
        print(f"N'={a.Np} {row['model']} seeds {s1},{s2}: a={aa:.4f} E[nu]={m[0]:.4f} E[nu^2]={m[1]:.4f} "
              f"(free {2 * aa * aa - aa ** 3:.4f}) E[nu^3]={m[2]:.4f} (free {5 * aa ** 3 - 6 * aa ** 4 + 2 * aa ** 5:.4f}) "
              f"KS_W(a,a)={row['ks_wachter_aa']:.3f}", flush=True)


if __name__ == "__main__":
    main()
