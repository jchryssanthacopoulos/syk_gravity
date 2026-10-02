#!/usr/bin/env python3
r"""
Exact size content of the chord-1 sector (first Wiener chaos) of the BPS projector, any (N', F, p).

The (1,1) chaos of the couplings carries Lambda^p (x) conj(Lambda^p) = (+)_{r=0..p} V_r, each irrep once, so the
chaos-1 component of P is P_[1] = sum_r g_r Y^(r) with Y = Q^dag Q restricted to Lambda^F (its size-r part, r >= 1,
is pure chaos 1; at half filling the even-r directions are H^(r) = 2 Y^(r), and odd r vanish).  Hence, exactly,
        J_1r = (E <P, Y^(r)>)^2 / (d E ||Y^(r)||^2)        (r >= 1),
estimated from per-realization overlaps.  We also report the size shares of H = {Q, Q^dag} (traceless part).
Uses: (i) odd N' (no particle-hole symmetry): is the top size r = p attained?  (ii) p = 5: does the chord-1 sector
have H's size proportions?  (iii) half filling, p = 3: control (J_11 = J_13 = 0, J_12 = pi_1 closed form).
Research note: research/notes/chord_size_dictionary_2026-10-02.md, Secs. 3-4.

Usage: .venv/bin/python scripts/run_chord1_sizes.py --Np 9 --F 4 --p 3 --samples 2000 --seed 5 \
           --out results/data/chord_chaos_2026-10-01/checks/chord1_p3_N9_F4.json
Memory: a few dense D x D matrices (N' <= 12: < 0.5 GB).  Even N' at half filling uses the even-size projector.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from chord_chaos import EvenSizeProjector, Sector, bps_projector_hodge, random_couplings  # noqa: E402
from size_decomposition import HopTable, size_components  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Np", type=int, required=True)
    ap.add_argument("--F", type=int, default=None, help="charge sector (default N'//2)")
    ap.add_argument("--p", type=int, default=3)
    ap.add_argument("--samples", type=int, default=500)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    F = a.Np // 2 if a.F is None else a.F
    S = Sector(a.Np, a.p, F=F)
    even = (2 * F == a.Np)
    if even:
        proj = EvenSizeProjector(a.Np, F)
        sizes = proj.ks
    else:
        T = HopTable(a.Np, F)
        sizes = list(range(min(F, a.Np - F) + 1))
    rng = np.random.default_rng(a.seed)
    tY, nY, nH = [], [], []
    t0 = time.time()
    for s in range(a.samples):
        c = random_couplings(rng, S.M)
        Qu = S.up.matrix(c)
        Y = (Qu.conj().T @ Qu).toarray()
        H = S.hamiltonian(c)
        P, _, d = bps_projector_hodge(S, c)
        if even:
            # half filling: odd sizes of P vanish; for even r the chaos-1 direction is H^(r) = 2 (Q^dag Q)^(r).
            # (Q^dag Q itself has odd sizes, so the even-size projector must be applied to H, not to Y.)
            cH, _ = proj.components(H)
            cY = cH
        else:
            cY, _ = size_components(Y, T)
            cH, _ = size_components(H, T)
        tY.append([np.vdot(P, x).real for x in cY])
        nY.append([np.vdot(x, x).real for x in cY])
        nH.append([np.vdot(x, x).real for x in cH])
    tY, nY, nH = np.array(tY), np.array(nY), np.array(nH)
    m_t, se_t, m_n = tY.mean(0), tY.std(0, ddof=1) / np.sqrt(len(tY)), nY.mean(0)
    hshare = nH.mean(0)[1:] / nH.mean(0)[1:].sum()
    rows = []
    for i, r in enumerate(sizes):
        if r == 0 or m_n[i] <= 1e-9 * m_n.sum():
            continue
        h = len(tY) // 2
        J = m_t[i] ** 2 / (d * m_n[i])
        Jsplit = tY[:h, i].mean() * tY[h:, i].mean() / (d * m_n[i])       # unbiased for (E t)^2
        rows.append(dict(r=int(r), J1r=float(J), J1r_split=float(Jsplit), overlap_mean=float(m_t[i]),
                         overlap_se=float(se_t[i]), sigma=float(abs(m_t[i]) / se_t[i])))
    tot = sum(x["J1r"] for x in rows)
    for x in rows:
        x["share"] = x["J1r"] / tot
    out = dict(Np=a.Np, F=F, p=a.p, samples=a.samples, seed=a.seed, d=int(d), a=d / S.D, pi1=tot, sectors=rows,
               H_traceless_size_shares={int(r): float(v) for r, v in zip(sizes[1:], hshare)}, secs=time.time() - t0)
    print(json.dumps(out, indent=1))
    if a.out:
        os.makedirs(os.path.dirname(a.out), exist_ok=True)
        with open(a.out, "w") as fh:
            json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
