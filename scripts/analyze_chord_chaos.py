#!/usr/bin/env python3
r"""
Analysis of scripts/run_chord_chaos.py output: exact chord-number distribution pi_n (Wiener chaos of P), its joint
distribution J_nk with U(N) size, the Krylov distribution kappa_n, and BLY's P_n (docs/derivations.md D4).

Fits.  With t = cos^2(m pi / K) (lags m and K - m pooled):
   G(t)      = sum_{n>=0} pi_n t^n          (P overlaps;  pi_0 = a, pi_1 closed form are checks)
   Gt(t)     = sum_{n>=1} pit_n t^n         (Ptil = chaos >= 2 part;  pit_1 = 0 is the test of P_[1] = gamma H^(2))
   Gt_k(t)   = sum_{n>=1} Jt_nk t^n         (size-k part of Ptil;  selection rule Jt_nk = 0 for k > (p-1) n)
Least squares on the sample mean, errors by bootstrap over samples.

Krylov.  kappa_n = a p_n(0)^2, p_n the orthonormal polynomials of the disorder-averaged spectral measure of H on
Lambda^F (pooled eigenvalues of all angles and samples, Stieltjes procedure).  kappa_0 = a; kappa_1 = pi_1 up to 1/M.

Usage: .venv/bin/python scripts/analyze_chord_chaos.py --Np 12 [--nmax 6] [--boot 400]
"""
import argparse
import glob
import json
import os
import sys
from math import comb, exp

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from dssyk_n2 import hh_length_distribution  # noqa: E402
from size_decomposition import character_coeffs, dim_Vk  # noqa: E402


def load(path):
    rows = [json.loads(l) for l in open(path)]
    K = rows[0]["K"]
    G = np.array([r["G"] for r in rows]); Gt = np.array([r["Gt"] for r in rows])
    Ct = np.array([r["Ct"] for r in rows]); counts = np.array(rows[0]["counts"])
    Cf = np.array([r["Cf"] for r in rows]) if rows[0].get("Cf") is not None else None
    Sk = np.array([r["Sk"] for r in rows]) if rows[0].get("Sk") is not None else None
    return rows, K, G, Gt, Ct, Cf, counts, Sk


def pool_lags(X, K):
    """Average lags m and K - m; returns array over m = 0..K//2 (leading axis = samples)."""
    out = []
    for m in range(K // 2 + 1):
        if m == 0 or 2 * m == K:
            out.append(X[:, m])
        else:
            out.append(0.5 * (X[:, m] + X[:, K - m]))
    return np.stack(out, axis=1)


def classes_to_sizes(C, counts, N, F):
    """Class sums (..., r) -> size overlaps (..., k) = c_k dim V_k with c from character_coeffs on class means."""
    sh = C.shape
    flat = C.reshape(-1, sh[-1])
    out = np.zeros_like(flat)
    dims = np.array([dim_Vk(N, k) for k in range(sh[-1])], dtype=float)
    for i, row in enumerate(flat):
        out[i] = character_coeffs(row / counts, N, F) * dims
    return out.reshape(sh)


def lsq(t, y, powers):
    A = np.stack([t ** n for n in powers], axis=1)
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return coef


def fit_with_boot(t, Y, powers, boot, rng):
    """Y: samples x len(t). Returns (coef, err, resid_chi) using the sample mean and bootstrap."""
    mean = Y.mean(0)
    coef = lsq(t, mean, powers)
    bs = []
    S = Y.shape[0]
    for _ in range(boot):
        idx = rng.integers(0, S, S)
        bs.append(lsq(t, Y[idx].mean(0), powers))
    bs = np.array(bs)
    A = np.stack([t ** n for n in powers], axis=1)
    sem = Y.std(0, ddof=1) / np.sqrt(S)
    ok = sem > 1e-12 * max(np.abs(mean).max(), 1e-300)
    chi = np.sqrt(np.mean(((A @ coef - mean)[ok] / sem[ok]) ** 2)) if ok.any() else 0.0
    return coef, bs.std(0, ddof=1), chi


def gls_fit(t, Y, powers, offset=None, ridge=1e-12):
    """GLS fit of the sample mean of Y (samples x len(t)) to sum_{n in powers} c_n t^n (+ known offset(t)), using the
    sample covariance of the mean.  Points with zero variance are dropped.  Returns (coef, err, chi2/dof)."""
    S = Y.shape[0]
    y = Y.mean(0) - (offset if offset is not None else 0.0)
    C = np.cov(Y.T) / S
    keep = np.diag(C) > 1e-30
    tt, y, C = t[keep], y[keep], C[np.ix_(keep, keep)]
    C = C + ridge * np.trace(C) / len(C) * np.eye(len(C))
    A = np.stack([tt ** n for n in powers], axis=1)
    W = np.linalg.inv(C)
    F = A.T @ W @ A
    coef = np.linalg.solve(F, A.T @ W @ y)
    err = np.sqrt(np.diag(np.linalg.inv(F)))
    r = y - A @ coef
    dof = max(len(tt) - len(powers), 1)
    return coef, err, float(r @ W @ r) / dof


def krylov_kappa(spec_files, a, nmax):
    """kappa_n = a p_n(0)^2 for the pooled spectral measure (uniform weights)."""
    ev = np.concatenate([np.load(f).ravel() for f in spec_files])
    ev = np.where(np.abs(ev) < 1e-9 * np.abs(ev).max(), 0.0, ev)
    scale = ev.mean()
    x = ev / scale
    w = np.full(len(x), 1.0 / len(x))
    # Stieltjes procedure
    p_prev = np.zeros_like(x); p = np.ones_like(x)
    nrm = np.sqrt(np.sum(w * p * p)); p /= nrm
    vals = [1.0 / nrm]           # p_0(0)
    b_prev = 0.0
    p0_prev, p0 = 0.0, 1.0 / nrm
    for n in range(nmax):
        alpha = np.sum(w * x * p * p)
        q = (x - alpha) * p - b_prev * p_prev
        b = np.sqrt(np.sum(w * q * q))
        q /= b
        p0_new = ((0.0 - alpha) * p0 - b_prev * p0_prev) / b
        p_prev, p, b_prev = p, q, b
        p0_prev, p0 = p0, p0_new
        vals.append(p0)
    vals = np.array(vals)
    atom = np.mean(ev == 0.0)
    return a * vals ** 2, atom, len(ev)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Np", type=int, required=True)
    ap.add_argument("--p", type=int, default=3)
    ap.add_argument("--nmax", type=int, default=6)
    ap.add_argument("--boot", type=int, default=400)
    ap.add_argument("--dir", default=os.path.join(ROOT, "results/data/chord_chaos_2026-10-01"))
    ap.add_argument("--json", default=None, help="write the fitted tables to this json file")
    a = ap.parse_args()
    N, p = a.Np, a.p
    F = N // 2
    rows, K, G, Gt, Ct, Cf, counts, Sk = load(os.path.join(a.dir, f"chord_chaos_p{p}_N{N}.jsonl"))
    S = len(rows)
    aa, pi1 = rows[0]["a"], rows[0]["pi1_closed"]
    rng = np.random.default_rng(12345)
    t = np.cos(np.pi * np.arange(K // 2 + 1) / K) ** 2
    Gp, Gtp = pool_lags(G, K), pool_lags(Gt, K)
    print(f"N'={N} p={p}: {S} samples x K={K} angles; a={aa:.6f}, closed-form pi_1={pi1:.6f}")
    print(f"  mean ||Ptil||^2/d = {Gtp[:,0].mean():.6f} +- {Gtp[:,0].std(ddof=1)/np.sqrt(S):.6f}"
          f"   (exact 1 - a - pi_1 = {1 - aa - pi1:.6f})")

    # --- total chord distribution from P overlaps (checks pi_0, pi_1) and from Ptil (pi_n, n >= 2)
    powsG = list(range(0, a.nmax + 1))
    cG, eG, chiG = fit_with_boot(t, Gp, powsG, a.boot, rng)
    powsT = list(range(1, a.nmax + 1))
    cT, eT, chiT = fit_with_boot(t, Gtp, powsT, a.boot, rng)
    print("  fit of G  (free pi_0..):  " + "  ".join(f"pi{n}={c:.4f}({e:.4f})" for n, c, e in zip(powsG, cG, eG))
          + f"   rms resid/sem={chiG:.2f}")
    print("  fit of Gt (free n>=1):    " + "  ".join(f"pt{n}={c:.4f}({e:.4f})" for n, c, e in zip(powsT, cT, eT))
          + f"   rms resid/sem={chiT:.2f}")
    # constrained: pit_1 = 0
    powsT2 = list(range(2, a.nmax + 1))
    cT2, eT2, chiT2 = fit_with_boot(t, Gtp, powsT2, a.boot, rng)
    print("  fit of Gt (n>=2 only):    " + "  ".join(f"pi{n}={c:.4f}({e:.4f})" for n, c, e in zip(powsT2, cT2, eT2))
          + f"   rms resid/sem={chiT2:.2f}")

    # --- preferred estimator: full-P overlaps with the exact n <= 1 part subtracted (P is scale invariant, so its
    #     overlaps fluctuate less than those of Ptil, whose gamma H^(2) term scales with |omega|^2)
    cR, eR, chiR = fit_with_boot(t, Gp - aa - pi1 * t[None, :], powsT2, a.boot, rng)
    print("  fit of G - a - pi1 t:     " + "  ".join(f"pi{n}={c:.4f}({e:.4f})" for n, c, e in zip(powsT2, cR, eR))
          + f"   rms resid/sem={chiR:.2f}   sum={cR.sum():.4f} (exact {1 - aa - pi1:.4f})")

    # --- size-resolved
    use_full = Cf is not None
    if Sk is not None:
        Jt = Sk.copy()                                 # exact Casimir size overlaps of the full P
        Jt[:, :, 0] -= aa
        Jt[:, :, 2] -= pi1 * (np.cos(np.pi * np.arange(K) / K) ** 2)[None, :]
        print("  size-resolved fits use EXACT (Casimir) size overlaps of P with the n <= 1 parts subtracted")
    elif use_full:
        Jt = classes_to_sizes(Cf, counts, N, F)       # full P; subtract exact n <= 1 parts below
        Jt[:, :, 0] -= aa
        Jt[:, :, 2] -= pi1 * (np.cos(np.pi * np.arange(K) / K) ** 2)[None, :]
        print("  size-resolved fits use FULL-P classes with J_0k = a delta_k0, J_1k = pi1 delta_k2 subtracted")
    else:
        Jt = classes_to_sizes(Ct, counts, N, F)       # Ptil classes
    rmax = Jt.shape[-1] - 1
    Jtp = np.stack([pool_lags(Jt[:, :, k], K) for k in range(rmax + 1)], axis=2)   # samples x m x k
    print(f"  chaos>=2 size weights at t=1 (k=0..{rmax}): " + " ".join(f"{x:.5f}" for x in Jtp[:, 0, :].mean(0)))
    Jfree = np.zeros((rmax + 1, len(powsT))); Jferr = np.zeros_like(Jfree)
    Jcon = np.zeros((rmax + 1, a.nmax + 1)); Jcerr = np.zeros_like(Jcon)
    for k in range(rmax + 1):
        c, e, _ = fit_with_boot(t, Jtp[:, :, k], powsT, a.boot, rng)
        Jfree[k], Jferr[k] = c, e
        nmin = max(2, -(-k // (p - 1)))           # selection rule: n >= k/(p-1)
        pw = list(range(nmin, a.nmax + 1))
        if k >= 1 and pw:
            c2, e2, _ = fit_with_boot(t, Jtp[:, :, k], pw, a.boot, rng)
            Jcon[k, nmin:], Jcerr[k, nmin:] = c2, e2
    if Sk is not None:
        # free fits (powers 0..nmax) of the exact size overlaps of the FULL P; n <= 1 coefficients are checks
        Skp = np.stack([pool_lags(Sk[:, :, k], K) for k in range(rmax + 1)], axis=2)
        Jfull = np.zeros((rmax + 1, a.nmax + 1)); Jfullerr = np.zeros_like(Jfull)
        for k in range(0, rmax + 1, 2):
            c, e, _ = fit_with_boot(t, Skp[:, :, k], list(range(0, a.nmax + 1)), a.boot, rng)
            Jfull[k], Jfullerr[k] = c, e
        print("  FREE fits of exact size overlaps (rows k even, cols n=0..nmax); checks: J_00=a, J_12=pi1:")
        for k in range(0, rmax + 1, 2):
            print(f"    k={k}: " + "  ".join(f"{c:+.4f}({e:.4f})" for c, e in zip(Jfull[k], Jfullerr[k])))
        pin_free = Jfull.sum(0)
        print("  pi_n from free size fits (sum over k): " + " ".join(f"{x:.4f}" for x in pin_free))
    if Sk is not None:
        # constrained GLS: exact J_0k = a delta_k0, J_1k = pi1 delta_k2 imposed; selection rule J_nk = 0 for n < k/(p-1)
        Skp2 = np.stack([pool_lags(Sk[:, :, k], K) for k in range(rmax + 1)], axis=2)
        Jg = np.zeros((rmax + 1, a.nmax + 1)); Jgerr = np.zeros_like(Jg)
        Jg[0, 0] = aa
        if rmax >= 2:
            Jg[2, 1] = pi1
        print("  CONSTRAINED GLS fits (exact n<=1 imposed; selection rule n >= k/(p-1)):")
        for k in range(2, rmax + 1, 2):
            nmin = max(2, -(-k // (p - 1)))
            pw = list(range(nmin, a.nmax + 1))
            off = pi1 * t if k == 2 else None
            c, e, chi = gls_fit(t, Skp2[:, :, k], pw, off)
            Jg[k, nmin:], Jgerr[k, nmin:] = c, e
            print(f"    k={k}: " + "  ".join(f"n{n}:{cc:+.4f}({ee:.4f})" for n, cc, ee in zip(pw, c, e)) + f"   chi2/dof={chi:.2f}")
        print("  pi_n (constrained GLS): " + " ".join(f"{x:.4f}" for x in Jg.sum(0)) +
              f"   shares J_n,2n/pi_n: " + " ".join(f"{Jg[2*n, n]/Jg[:, n].sum():.3f}" for n in range(2, min(a.nmax, rmax // 2) + 1)
                                                  if abs(Jg[:, n].sum()) > 1e-9))
    print("  Jt_nk, unconstrained fit (rows k, cols n=1..nmax):")
    for k in range(rmax + 1):
        print(f"    k={k}: " + "  ".join(f"{c:+.4f}({e:.4f})" for c, e in zip(Jfree[k], Jferr[k])))
    # assemble the joint distribution J_nk (n=0,1 exact)
    J = Jcon.copy(); J[0, 0] = aa; J[2, 1] = pi1
    Jerr = Jcerr.copy()
    print("  Joint J_nk (rows k, cols n=0..nmax), n<=1 exact, n>=2 fitted with selection rule k <= (p-1) n:")
    for k in range(rmax + 1):
        if k % 2 and np.all(np.abs(J[k]) < 5e-4):
            continue
        print(f"    k={k}: " + "  ".join(f"{c:.4f}({e:.4f})" if c != 0 else "   0         "
                                        for c, e in zip(J[k], Jerr[k])))
    pin = J.sum(0); pinerr = np.sqrt((Jerr ** 2).sum(0))
    wk = J.sum(1)
    lam = 2 * p * p / N
    Pn = hh_length_distribution(lam, 0.0)
    specs = sorted(glob.glob(os.path.join(a.dir, f"spectra_p{p}_N{N}_s*.npy")))
    kap, atom, nev = krylov_kappa(specs, aa, a.nmax)
    print(f"  Krylov: {nev} pooled eigenvalues, atom at 0 = {atom:.6f} (a = {aa:.6f})")
    print("   n   pi_n (exact chaos)    kappa_n (Krylov)   P_n (BLY, lam=2p^2/N)   w_{(p-1)n} (size)")
    for n in range(a.nmax + 1):
        k = (p - 1) * n
        print(f"   {n}   {pin[n]:.4f} ({pinerr[n]:.4f})     {kap[n]:.4f}            {Pn[n]:.4f}"
              f"                {wk[k] if k <= rmax else float('nan'):.4f}")
    if Sk is not None:
        # compact dictionary table from the FREE exact-size fits (n <= nrep reported; higher n are cutoff-sensitive)
        from size_decomposition import chi_hat
        chi = np.array([chi_hat(N, 1, k) for k in range(rmax + 1)])
        x_ds = exp(-2 * p / N)
        nrep = min(4, a.nmax)
        print("  DICTIONARY (free exact-size fits):")
        print("   n  pi_n     kappa_n  P_n(BLY)  J_n,(p-1)n  share   w_(p-1)n  chibar_n  x^2n(DS)  (1-2p/N)^2n")
        Ttot = 0.0
        for n in range(0, a.nmax + 1):
            pin_ = Jfull[:, n].sum()
            k = (p - 1) * n
            Jd = Jfull[k, n] if k <= rmax else float('nan')
            chib = (Jfull[:, n] @ chi) / pin_ if abs(pin_) > 1e-12 else float('nan')
            Ttot += Jfull[:, n] @ chi
            if n <= nrep:
                print(f"   {n}  {pin_:.4f}   {kap[n]:.4f}   {Pn[n]:.4f}    {Jd:.4f}      {Jd/pin_ if pin_ else float('nan'):.3f}   "
                      f"{(Jfull[k].sum() if k <= rmax else float('nan')):.4f}    {chib:+.4f}   {x_ds**(2*n):.4f}    {(1-2*p/N)**(2*n):.4f}")
        mean_n = sum(n * Jfull[:, n].sum() for n in range(a.nmax + 1))
        mean_k = sum(k * Jfull[k].sum() for k in range(rmax + 1))
        print(f"   <n>_pi = {mean_n:.4f}   <k>_w = {mean_k:.4f}   <k>/<n> = {mean_k/mean_n:.3f} (bound p-1 = {p-1});"
              f"  T2 = sum J chi = {Ttot:.4f}   BLY sum P_n x^2n = {sum(Pn[n]*x_ds**(2*n) for n in range(len(Pn))):.4f}")
    if a.json:
        out = dict(Np=N, p=p, samples=S, K=K, a=aa, pi1_closed=pi1, J=J.tolist(), Jerr=Jerr.tolist(),
                   J_free=(Jfull.tolist() if Sk is not None else None),
                   J_gls=(Jg.tolist() if Sk is not None else None), J_gls_err=(Jgerr.tolist() if Sk is not None else None),
                   J_free_err=(Jfullerr.tolist() if Sk is not None else None),
                   J_unconstrained=Jfree.tolist(), J_unconstrained_err=Jferr.tolist(), pi=pin.tolist(),
                   pi_err=pinerr.tolist(), kappa=kap.tolist(), P_BLY=Pn[:a.nmax + 1].tolist(), w_from_J=wk.tolist(),
                   fitG=cG.tolist(), fitG_err=eG.tolist(), fitGt=cT.tolist(), fitGt_err=eT.tolist())
        with open(a.json, "w") as fh:
            json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
