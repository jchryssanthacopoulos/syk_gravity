#!/usr/bin/env python3
r"""
q_scan_mf.py -- MATRIX-FREE estimate of the decoder freeness failure delta m_2,
for (N,p,q) too large to store the BPS basis.

Idea: never build P_B or its basis. The harmonic projector is P_B = projector onto
ker L, L = Q^dag Q + Q Q^dag the wedge-Laplacian of the q-form supercharge Q = C wedge,
a C x C SPARSE Hermitian operator (C = binom(N',p)). We apply P_B to a vector by one
iterative solve:
    (I - P_B) v = argmin_x || L x - L v ||   (min-norm, via CG from x0=0),   P_B v = v - x.
Then Hutchinson sampling with random complex z (E[z z^dag] = I) gives the traces
    d          = Tr P_B                  = E[ z^dag P_B z ],
    Tr(P_B Pi) = E[ z^dag P_B Pi_T z ],
    Tr((P_B Pi)^2) = E[ z^dag P_B Pi_T P_B Pi_T z ],
hence  m_k^meas = Tr((P_B Pi)^k)/d  and  delta m_k = m_k^meas - m_k^Wachter(a,b).
Memory is O(C) per vector (a handful of vectors + the sparse L), so this reaches
N ~ 20-24 where basis storage (O(C*d)) is hopeless.

Trade-off: stochastic noise. Error on delta m_2 falls like 1/sqrt(K*batches); the run
prints a jackknife error bar so you know when it is tight enough (aim < 0.003 to read
the trend). Validate against exact q_scan.py at N=8,10 before trusting large N.

USAGE
    python q_scan_mf.py --N 10 --p 5 --q 3 --K 400            # single point + error bar
    python q_scan_mf.py --N 10 --p 5 --q 3 --K 400 --check    # compare to exact basis
    python q_scan_mf.py --N 16 --p 8 --q 4 --K 600            # large N, matrix-free
Requires: numpy, scipy only.
"""
import argparse, inspect
from math import comb
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import cg
from q_scan import subset_index, form_wedge, rand_form, wachter_moment, harmonic

_CGKW = "rtol" if "rtol" in inspect.signature(cg).parameters else "tol"


def build_L(C, n, p, q):
    """wedge-Laplacian on Lambda^p: L = Wup^dag Wup + Wdn Wdn^dag, ker L = harmonic space."""
    Wup = form_wedge(C, n, p, q)                       # Lambda^p -> Lambda^{p+q}
    L = (Wup.conj().T @ Wup).tocsr()
    if p - q >= 0:
        Wdn = form_wedge(C, n, p - q, q)               # Lambda^{p-q} -> Lambda^p
        L = (L + (Wdn @ Wdn.conj().T)).tocsr()
    return L


def real_embed(L):
    """real 2C x 2C symmetric embedding of a complex Hermitian L (same spectrum, x2)."""
    A, B = L.real, L.imag
    return sp.bmat([[A, -B], [B, A]], format="csr")


def make_PB(L, tol, maxiter):
    """return applier v |-> P_B v = v - (min-norm solution of L x = L v)."""
    Lr = real_embed(L)
    Cn = L.shape[0]

    def PB(v):
        b = L @ v
        br = np.concatenate([b.real, b.imag])
        xr, _ = cg(Lr, br, maxiter=maxiter, **{_CGKW: tol})
        return v - (xr[:Cn] + 1j * xr[Cn:])
    return PB


def estimate(N, p, q, seed, K, batches, cg_tol, cg_iter):
    n, D = N + 1, comb(N + 1, p)
    masks, _ = subset_index(n, p)
    slot = np.array([0.0 if (mm >> N) & 1 else 1.0 for mm in masks])
    b = slot.sum() / D
    L = build_L(rand_form(seed, n, q), n, p, q)
    PB = make_PB(L, cg_tol, cg_iter)
    rng = np.random.default_rng(1000 + seed)

    # accumulate per-batch trace sums; pool them for the (low-bias) central value and
    # jackknife over batches for the error bar. Pooling makes the ratio-estimator bias
    # O(1/[K*batches]) instead of O(1/K).
    SPB, S1, S2 = [], [], []
    for _ in range(batches):
        s_pb = s_1 = s_2 = 0.0
        for _ in range(K):
            z = (rng.standard_normal(D) + 1j * rng.standard_normal(D)) / np.sqrt(2)
            pz = PB(z)                        # P_B z
            s_pb += np.vdot(z, pz).real       # -> Tr P_B = d
            py = PB(slot * z)                 # P_B Pi_T z
            s_1 += np.vdot(z, py).real        # -> Tr(P_B Pi_T)
            p4 = PB(slot * py)                # P_B Pi_T P_B Pi_T z
            s_2 += np.vdot(z, p4).real        # -> Tr((P_B Pi_T)^2)
        SPB.append(s_pb); S1.append(s_1); S2.append(s_2)
    SPB, S1, S2 = np.array(SPB), np.array(S1), np.array(S2)

    def deltas(spb, s1, s2, cnt):
        d = spb / cnt; a = d / D
        return (a, s1 / spb - wachter_moment(a, b, 1),
                s2 / spb - wachter_moment(a, b, 2))

    tot = K * batches
    a, dm1, dm2 = deltas(SPB.sum(), S1.sum(), S2.sum(), tot)   # pooled central value
    j1, j2 = [], []                                            # leave-one-batch-out
    for j in range(batches):
        keep = np.arange(batches) != j
        _, d1, d2 = deltas(SPB[keep].sum(), S1[keep].sum(), S2[keep].sum(),
                           (batches - 1) * K)
        j1.append(d1); j2.append(d2)
    j1, j2 = np.array(j1), np.array(j2)
    jse = lambda jk: np.sqrt((batches - 1) / batches * np.sum((jk - jk.mean()) ** 2)) \
        if batches > 1 else float("nan")
    return dict(a=a, b=b, D=D, dm1=dm1, dm1_se=jse(j1), dm2=dm2, dm2_se=jse(j2))


def exact_point(N, p, q, seed):
    """exact basis result for validation (small N only)."""
    n, D = N + 1, comb(N + 1, p)
    masks, _ = subset_index(n, p)
    slot = np.array([0.0 if (mm >> N) & 1 else 1.0 for mm in masks])
    b = slot.sum() / D
    B = harmonic(rand_form(seed, n, q), n, p, q)
    d = B.shape[1]
    a = d / D
    KT = B.conj().T @ (slot[:, None] * B)
    lam = np.clip(np.linalg.eigvalsh((KT + KT.conj().T) / 2).real, 0, 1)
    return a, b, (lam ** 2).mean() - wachter_moment(a, b, 2)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--p", type=int, required=True)
    ap.add_argument("--q", type=int, required=True)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--K", type=int, default=300, help="Hutchinson samples per batch")
    ap.add_argument("--batches", type=int, default=8, help="batches (for the error bar)")
    ap.add_argument("--cg-tol", type=float, default=1e-9)
    ap.add_argument("--cg-iter", type=int, default=5000)
    ap.add_argument("--check", action="store_true", help="also compute exact (small N)")
    a = ap.parse_args()

    r = estimate(a.N, a.p, a.q, a.seed, a.K, a.batches, a.cg_tol, a.cg_iter)
    print(f"(N={a.N}->{a.N+1}, p={a.p}, q={a.q})  D={r['D']}  a~{r['a']:.3f}  b={r['b']:.3f}"
          f"   [{a.K*a.batches} samples]")
    print(f"  matrix-free:  delta m_2 = {r['dm2']:+.4f} +/- {r['dm2_se']:.4f}"
          f"     (delta m_1 = {r['dm1']:+.4f} +/- {r['dm1_se']:.4f})")
    if a.check:
        ea, eb, edm2 = exact_point(a.N, a.p, a.q, a.seed)
        z = (r["dm2"] - edm2) / (r["dm2_se"] + 1e-12)
        print(f"  exact basis:  delta m_2 = {edm2:+.4f}   (a={ea:.3f})   "
              f"[matrix-free off by {z:+.1f} sigma]")


if __name__ == "__main__":
    main()
