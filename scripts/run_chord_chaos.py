#!/usr/bin/env python3
r"""
Joint (chord number n, U(N) size k) distribution of the BPS projector of one-flavor N=2 SYK (docs/derivations.md D4).

Chord number = Wiener-chaos degree of P(omega) in the Gaussian couplings.  For OU-correlated couplings
omega_x = x omega + sqrt(1 - x^2) omega',   E tau(P(omega) P(omega_x)) = sum_n pi_n x^(2n)   and, size-resolved,
E <P^(k)(omega), P^(k)(omega_x)> / d = sum_n J_nk x^(2n).

Angle trick: per sample draw omega, omega'; for theta_j = j pi / K (j = 0..K-1) build the BPS projector P_j of
omega_j = cos(theta_j) omega + sin(theta_j) omega'.  Every pair (i, j) is an OU pair with x = cos((j - i) pi / K).
Variance reduction (p = 3, half filling): Ptil = P - a 1 - gamma (H - Tr H / D) is the chaos >= 2 part of P exactly
(P_[0] = a 1, P_[1] = gamma H^(2), src/chord_chaos.py), so E tau(Ptil Ptil_x) = sum_{n>=2} pi_n x^(2n).

Per sample and lag m = (j - i) mod K, averaged over the K pairs at that lag, the script records
  G[m]      = Tr(P_i P_j) / d
  Gt[m]     = Tr(Ptil_i Ptil_j) / d
  Ct[m][r]  = class sums over basis pairs (S, T) with |S \ T| = r of Re(conj(Ptil_i)_ST (Ptil_j)_ST) / d
  C[m][r]   = same for P (only if --full-classes)
and the spectra of H (for the Krylov comparison).  Analysis: scripts/analyze_chord_chaos.py.

Memory: K dense D x D complex matrices for Ptil (+ K for P with --full-classes) + eigh workspace.
  N'=12: ~0.5 GB (K=16);  N'=14: ~16 * 0.19 GB * (1 or 2) + 1 GB.  Use --dry-run for the estimate.
Run heavy jobs under scripts/memwatch.py, e.g.
  .venv/bin/python scripts/memwatch.py --limit-gb 4 --log results/data/chord_chaos_2026-10-01/memwatch.jsonl -- \
      .venv/bin/python scripts/run_chord_chaos.py --Np 12 --K 16 --samples 40 --seed0 0
"""
import argparse
import json
import os
import subprocess
import sys
import time
from math import comb

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from chord_chaos import EvenSizeProjector, Sector, bps_basis_from_H, bps_projector_hodge, random_couplings  # noqa: E402


def git_state():
    try:
        g = subprocess.run(["git", "-C", ROOT, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", ROOT, "status", "--porcelain"], capture_output=True, text=True).stdout.strip()
        return g + ("-dirty" if dirty else "")
    except OSError:
        return "unknown"


def mem_estimate_gb(N, K, full, exact=False, ptil=True):
    D = comb(N, N // 2)
    per = D * D * 16 / 1e9
    return (K * per * ((1 if ptil else 0) + (1 if full else 0)) + (K * per + 8 * per if exact else 0)
            + 6 * per + 0.3)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Np", type=int, required=True, help="number of complex modes N' (even: half filling)")
    ap.add_argument("--p", type=int, default=3, help="fermions per supercharge term")
    ap.add_argument("--K", type=int, default=16, help="number of rotation angles per sample")
    ap.add_argument("--samples", type=int, default=10)
    ap.add_argument("--seed0", type=int, default=0, help="sample s uses numpy seed 1_000_003 * (seed0 + s) + 17")
    ap.add_argument("--full-classes", action="store_true", help="also record size classes of P itself")
    ap.add_argument("--exact-sizes", action="store_true",
                    help="record exact (Casimir) size-resolved overlaps Sk[m][k] = mean_i Re<P_i^(k), P_j>/d")
    ap.add_argument("--no-torus", action="store_true", help="skip the torus class sums (Ct, Cf)")
    ap.add_argument("--no-ptil", action="store_true",
                    help="do not store Ptil = P - a - gamma H^(2) (skips Gt; halves memory at large N')")
    ap.add_argument("--out", default=os.path.join(ROOT, "results/data/chord_chaos_2026-10-01"))
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    est = mem_estimate_gb(a.Np, a.K, a.full_classes, a.exact_sizes, not a.no_ptil)
    print(f"N'={a.Np} p={a.p} K={a.K}: estimated peak memory {est:.2f} GB", flush=True)
    if a.dry_run:
        return
    os.makedirs(a.out, exist_ok=True)
    S = Sector(a.Np, a.p)
    d_exp = None
    started = time.strftime("%Y%m%dT%H%M%S")
    meta = dict(args=vars(a), git=git_state(), started=started, D=S.D, M=S.M, F=S.F,
                E_trH=S.E_trH(), E_normH_nonid=S.E_normH_nonid(), mem_estimate_gb=est,
                description="angle-rotated OU pairs; chord number = Wiener chaos degree; see docs/derivations.md D4")
    tag = f"p{a.p}_N{a.Np}"
    with open(os.path.join(a.out, f"meta_{tag}_{started}.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    thetas = np.pi * np.arange(a.K) / a.K
    I = np.eye(S.D)
    if not a.no_torus:
        _ = S.R  # build class index
    proj = EvenSizeProjector(a.Np, S.F) if a.exact_sizes else None
    for s in range(a.samples):
        t0 = time.time()
        seed = 1_000_003 * (a.seed0 + s) + 17
        rng = np.random.default_rng(seed)
        c1, c2 = random_couplings(rng, S.M), random_couplings(rng, S.M)
        Bs, Pt, Pf, spectra, gaps, Pd = [], [], [], [], [], ([] if proj is not None else None)
        for th in thetas:
            c = np.cos(th) * c1 + np.sin(th) * c2
            H = S.hamiltonian(c)
            if proj is not None:                       # Hodge projector: two thin SVDs (no full eigensolve)
                Pm, w, d = bps_projector_hodge(S, c, d_exp)
                gap = np.inf
            else:
                B, w, gap = bps_basis_from_H(H, d_exp)
                d = B.shape[1]
                Pm = B @ B.conj().T
                Bs.append(B)
            if d_exp is None:
                d_exp = d
                pi1, gamma = S.pi1_closed(d)
                aa = d / S.D
            trH = np.trace(H).real
            if not a.no_ptil:
                Pt.append(Pm - aa * I - gamma * (H - (trH / S.D) * I))
            if a.full_classes:
                Pf.append(Pm)
            if Pd is not None:
                Pd.append(Pm)
            spectra.append(w)
            gaps.append(gap)
            del H, Pm
        K = a.K
        G = np.zeros(K); Gt = np.zeros(K)
        Ct = np.zeros((K, S.rmax + 1)); Cf = np.zeros((K, S.rmax + 1))
        nks = len(proj.ks) if proj is not None else 0
        Sk = np.zeros((K, S.rmax + 1)); size_resid = 0.0
        for i in range(K):
            if proj is not None:
                comps, res = proj.components(Pd[i], check=(i == 0))
                if res is not None:
                    size_resid = res
            for m in range(K):
                j = (i + m) % K
                if Pd is not None:
                    G[m] += np.vdot(Pd[i], Pd[j]).real / d_exp
                else:
                    G[m] += np.linalg.norm(Bs[i].conj().T @ Bs[j]) ** 2 / d_exp
                if not a.no_ptil:
                    Gt[m] += np.vdot(Pt[i], Pt[j]).real / d_exp
                if not a.no_torus:
                    Ct[m] += S.class_sums(Pt[i], Pt[j]) / d_exp
                    if a.full_classes:
                        Cf[m] += S.class_sums(Pf[i], Pf[j]) / d_exp
                if proj is not None:
                    for ii, k in enumerate(proj.ks):
                        Sk[m, k] += np.vdot(comps[ii], Pd[j]).real / d_exp
            if proj is not None:
                del comps
        G /= K; Gt /= K; Ct /= K; Cf /= K; Sk /= K
        row = dict(Np=a.Np, p=a.p, K=K, sample=a.seed0 + s, seed=seed, d=d_exp, a=aa, pi1_closed=pi1, gamma=gamma,
                   G=G.tolist(), Gt=Gt.tolist(), Ct=Ct.tolist(), counts=S.class_counts().tolist(),
                   Cf=(Cf.tolist() if (a.full_classes and not a.no_torus) else None),
                   Sk=(Sk.tolist() if proj is not None else None), size_resid=size_resid,
                   min_gap_ratio=float(np.min(gaps)),
                   secs=time.time() - t0)
        with open(os.path.join(a.out, f"chord_chaos_{tag}.jsonl"), "a") as fh:
            fh.write(json.dumps(row) + "\n")
        np.save(os.path.join(a.out, f"spectra_{tag}_s{a.seed0 + s}.npy"), np.array(spectra, dtype=np.float64))
        print(f"sample {a.seed0 + s}: Gt[0]={Gt[0]:.5f} (exact mean {1 - aa - pi1:.5f})  G[K/2]={G[K // 2]:.5f} "
              f"(a={aa:.5f})  {row['secs']:.1f}s", flush=True)
        del Bs, Pt, Pf, Pd


if __name__ == "__main__":
    main()
