#!/usr/bin/env python3
r"""
q_scan.py -- does turning the supercharge degree q into a knob suppress the
decoder's freeness failure?

For the uplift decoder N -> N'=N+1 with a *q-form* supercharge, we build
  P_B  = harmonic p-forms of the enlarged theory  (rank d, ker(C wedge) cap ker(C wedge)^dag,
         C a generic complex q-form),
  Pi_T = slot projector (p-forms without the new mode),               (rank m = C(N,p))
and study O_T = P_B Pi_T P_B, whose d eigenvalues should approach the Wachter law
Bern(a) |x| Bern(b),  a = d/D,  b = m/D,  D = C(N',p).

The freeness failure is measured *per BPS state* by
  delta m_k = m_k^meas - m_k^Wachter,      m_k^meas = (1/d) Tr(O_T^k),
whose leading, dominant piece is t_2 ~ delta m_2 (see the draft). If a large-q or
double-scaling limit restores freeness, t_2 must fall with q. This script tabulates
t_2 (and delta m_k) across a grid of (N, p, q) so you can see the trend directly.

Everything is cheap: only the d-dimensional compression K = B^dag Pi_T B is
diagonalized (no D x D operators), and the Wachter moments are closed form.

USAGE
    python q_scan.py                         # default grid (fixed (N,p), vary q)
    python q_scan.py --N 8 10 12 --q 2 3 4 5 --seeds 2
    python q_scan.py --half                  # half-filling p=floor((N+1)/2) only
    python q_scan.py --N 10 --p 5 --q 3 4 5  # explicit

  *** the clean double-scaling test ***  (needs larger N; run on real hardware)
    python q_scan.py --target-a 0.5 --q 3 4 5 6 --Nmax 20
      For each q it searches (N,p) for the one whose BPS fraction a=d/D is closest
      to the target, and reports delta m_2 along that FIXED-a ray. If delta m_2 keeps
      falling with q at fixed a, the freeness restoration is genuine (not the trivial
      a->1, P_B->identity effect that contaminates the fixed-(N,p) scan).

Requires: numpy, scipy.  (pip install numpy scipy)
"""
import argparse, itertools
from math import comb
import numpy as np
import scipy.sparse as sp


# ----------------------------------------------------------------- combinatorics
def subset_index(n, p):
    if p < 0 or p > n:
        return [], {}
    masks = []
    for S in itertools.combinations(range(n), p):
        m = 0
        for i in S:
            m |= (1 << i)
        masks.append(m)
    return masks, {m: i for i, m in enumerate(masks)}


def merge_sign(Tmask, Smask):
    """sign of e_T ^ e_S -> e_{sorted(T u S)} for disjoint bitmasks T, S."""
    s, t = 0, Tmask
    while t:
        b = t & (-t)
        s += bin(Smask & (b - 1)).count("1")
        t ^= b
    return -1.0 if (s & 1) else 1.0


def form_wedge(coeffs, n, p, deg):
    """sparse matrix of (omega ^ .): Lambda^p -> Lambda^{p+deg}, omega={tuple:val}."""
    rows = comb(n, p + deg) if 0 <= p + deg <= n else 0
    colsN = comb(n, p) if 0 <= p <= n else 0
    if rows == 0 or colsN == 0:
        return sp.csr_matrix((rows, colsN), dtype=complex)
    cols, _ = subset_index(n, p)
    _, rmap = subset_index(n, p + deg)
    I, J, V = [], [], []
    for cj, Smask in enumerate(cols):
        for T, c in coeffs.items():
            Tmask, ok = 0, True
            for i in T:
                if Smask & (1 << i):
                    ok = False
                    break
                Tmask |= (1 << i)
            if ok:
                I.append(rmap[Smask | Tmask]); J.append(cj)
                V.append(c * merge_sign(Tmask, Smask))
    return sp.csr_matrix((V, (I, J)), shape=(rows, colsN), dtype=complex)


def rand_form(seed, n, deg):
    rng = np.random.default_rng(seed)
    return {T: rng.standard_normal() + 1j * rng.standard_normal()
            for T in itertools.combinations(range(n), deg)}


def nullspace_dense(M):
    A = M.toarray()
    if min(A.shape) == 0:
        return np.eye(A.shape[1], dtype=complex)
    _, s, Vh = np.linalg.svd(A, full_matrices=True)
    tol = max(A.shape) * (s[0] if len(s) else 1.0) * np.finfo(float).eps
    rank = int((s > tol).sum())
    return Vh[rank:].conj().T


def harmonic(C, n, p, q):
    """orthonormal basis of B_n^p = ker(C^) cap ker(C^)^dag with a q-form C."""
    if p < 0 or p > n:
        return np.zeros((max(comb(n, p) if 0 <= p <= n else 0, 0), 0), dtype=complex)
    blocks = [form_wedge(C, n, p, q)]
    if p - q >= 0:
        blocks.append(form_wedge(C, n, p - q, q).conj().transpose())
    B = nullspace_dense(sp.vstack(blocks).tocsr())
    if B.shape[1]:
        B, _ = np.linalg.qr(B)
    return B


# ----------------------------------------------------------------- Wachter moments
def manova_moment(a, b, k):
    """m_k^MANOVA (D-normalized) = (1/k)[z^{k-1}] ((z+a)(z+b)/(1+z))^k."""
    num = np.array([a * b, a + b, 1.0])
    inv = np.array([(-1.0) ** j for j in range(k + 1)])
    h = np.convolve(num, inv)[:k + 1]
    p = np.array([1.0])
    for _ in range(k):
        p = np.convolve(p, h)[:k + 1]
    return p[k - 1] / k


def wachter_moment(a, b, k):
    return manova_moment(a, b, k) / a


# ----------------------------------------------------------------- one (N,p,q)
def measure(N, p, q, kmax, seeds):
    Np, D = N + 1, comb(N + 1, p)
    masks, _ = subset_index(Np, p)
    slot = np.array([0.0 if (mm >> N) & 1 else 1.0 for mm in masks])
    m = int(slot.sum()); b = m / D
    dms, d_last, a_last = [], 0, 0.0
    for s in seeds:
        C = rand_form(s, Np, q)
        B = harmonic(C, Np, p, q)
        d = B.shape[1]
        if d == 0:
            continue
        a = d / D
        KT = B.conj().T @ (slot[:, None] * B)
        lam = np.clip(np.linalg.eigvalsh((KT + KT.conj().T) / 2).real, 0.0, 1.0)
        mmeas = np.array([(lam ** k).mean() for k in range(1, kmax + 1)])
        mwach = np.array([wachter_moment(a, b, k) for k in range(1, kmax + 1)])
        dms.append(mmeas - mwach)
        d_last, a_last = d, a
    if not dms:
        return None
    dm = np.mean(dms, axis=0)
    return dict(N=N, p=p, q=q, D=D, d=d_last, a=a_last, b=b, dm=dm)


# ----------------------------------------------------------------- fixed-a ray
def fixed_a_ray(target_a, qs, Nmax, seeds, kmax, target_b=0.5):
    """For each q, find (N,p) whose (a,b) is closest to (target_a, target_b), and
    report delta m_k. Matching a ALONE is not enough: delta m depends on a, b, and N,
    so the trivial a->1 (P_B->identity) confound and the b-dependence must both be
    controlled. This routine matches (a,b) jointly and prints the residual mismatch
    so you can judge how clean each rung of the ray is. A decline of delta m_2 down
    the column at genuinely matched (a,b) is real large-q freeness restoration."""
    print(f"Double-scaling ray, target (a,b) = ({target_a}, {target_b}):")
    print(f"{'q':>3} {'N':>3} {'p':>3} {'D':>6} {'d':>6} {'a':>7} {'b':>7} "
          f"{'|da|+|db|':>9} {'dm2':>9} {'dm3':>9}")
    print("-" * 78)
    prev = None
    for q in qs:
        best, bestdist = None, 1e9
        for N in range(max(2 * q - 1, q + 1), Nmax + 1):
            h = (N + 1) // 2
            for p in range(max(q, h - 4), min(N + 1 - q, h + 4) + 1):
                r = measure(N, p, q, kmax, seeds)
                if r is None:
                    continue
                dist = abs(r["a"] - target_a) + abs(r["b"] - target_b)
                if dist < bestdist:
                    best, bestdist = r, dist
        if best is None:
            print(f"{q:>3}   (no BPS states found for N<= {Nmax})")
            continue
        dm = best["dm"]
        tag = ""
        if prev is not None and bestdist < 0.06:   # only compare well-matched rungs
            tag = " DOWN" if dm[1] < prev else (" UP" if dm[1] > prev else " flat")
        flag = "" if bestdist < 0.06 else "  (poor match, raise --Nmax)"
        print(f"{q:>3} {best['N']:>3} {best['p']:>3} {best['D']:>6} {best['d']:>6} "
              f"{best['a']:>7.3f} {best['b']:>7.3f} {bestdist:>9.3f} {dm[1]:>+9.4f} "
              f"{dm[2]:>+9.4f}{tag}{flag}")
        if bestdist < 0.06:
            prev = dm[1]
    print("\nRead: compare dm2 only across rungs with small |da|+|db| (well matched).")
    print("Falling dm2 at matched (a,b) = genuine large-q freeness restoration; flat or")
    print("rising = the fixed-(N,p) trend was the trivial a->1 (P_B->identity) effect.")
    print("Reaching matched (a,b) at q=5,6 needs N ~ 16-20, so raise --Nmax on real hardware.")


# ----------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--N", type=int, nargs="+", default=[8, 10, 12])
    ap.add_argument("--q", type=int, nargs="+", default=[2, 3, 4, 5])
    ap.add_argument("--p", type=int, nargs="+", default=None,
                    help="explicit p values; default scans p near half-filling")
    ap.add_argument("--half", action="store_true", help="only p=floor((N+1)/2)")
    ap.add_argument("--seeds", type=int, default=2)
    ap.add_argument("--kmax", type=int, default=6)
    ap.add_argument("--target-a", type=float, default=None,
                    help="double-scaling ray: for each q, pick (N,p) with (a,b) closest to target")
    ap.add_argument("--target-b", type=float, default=0.5,
                    help="target filling fraction b for the double-scaling ray (default 0.5)")
    ap.add_argument("--Nmax", type=int, default=14, help="max N searched in --target-a mode")
    a = ap.parse_args()

    if a.target_a is not None:
        fixed_a_ray(a.target_a, a.q, a.Nmax, range(a.seeds), a.kmax, a.target_b)
        return

    print(f"{'N':>3} {'p':>3} {'q':>3} {'D':>6} {'d':>6} {'a':>7} {'b':>7} "
          f"{'lambda=q^2/N':>12} {'t2~dm2':>9} {'dm3':>9} {'dm4':>9}")
    print("-" * 86)
    rows = []
    for N in a.N:
        if a.p is not None:
            ps = a.p
        elif a.half:
            ps = [(N + 1) // 2]
        else:
            h = (N + 1) // 2
            ps = [h - 1, h, h + 1]
        for p in ps:
            if not (0 <= p <= N + 1):
                continue
            for q in a.q:
                if q < 1 or q > N + 1:
                    continue
                r = measure(N, p, q, a.kmax, range(a.seeds))
                if r is None:
                    print(f"{N:>3} {p:>3} {q:>3} {comb(N+1,p):>6} {'0':>6} "
                          f"{'--':>7} {'--':>7} {q*q/N:>12.3f}   (no BPS states)")
                    continue
                rows.append(r)
                dm = r["dm"]
                print(f"{N:>3} {p:>3} {q:>3} {r['D']:>6} {r['d']:>6} "
                      f"{r['a']:>7.3f} {r['b']:>7.3f} {q*q/N:>12.3f} "
                      f"{dm[1]:>+9.4f} {dm[2]:>+9.4f} {dm[3]:>+9.4f}")
        print("-" * 86)

    # trend summary: does t2 fall with q at fixed (N,p)?
    print("\nTREND: t2 (~delta m_2) vs q at fixed (N,p)  [decreasing => toward freeness]")
    from collections import defaultdict
    byNp = defaultdict(dict)
    for r in rows:
        byNp[(r["N"], r["p"])][r["q"]] = r["dm"][1]
    for (N, p), qd in sorted(byNp.items()):
        qs = sorted(qd)
        if len(qs) < 2:
            continue
        seq = "  ".join(f"q{q}:{qd[q]:+.4f}" for q in qs)
        arrow = "DOWN" if qd[qs[-1]] < qd[qs[0]] else ("UP" if qd[qs[-1]] > qd[qs[0]] else "flat")
        print(f"  (N={N:2d}, p={p}):  {seq}   -> {arrow}")


if __name__ == "__main__":
    main()
