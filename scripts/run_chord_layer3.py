#!/usr/bin/env python3
r"""
Layer 3 of the resonant chord problem (docs/derivations.md D11 §4): the exact decoder of the resonant sector from the
first-order clouds plus the back-reaction operator R, and its free (Haar) chord model.

Exact structure (derived, D11 §4).
- Self-dressed h block (D6 §5): G = 1 + K^dag (1+R)^-1 K, N = 1 - K^dag (1-R)(1+R)^-2 K. In the edge-mode basis
  (singular basis of Q0: coexact(l'-1) -> exact(l'+2q-1)), K = e / s (first-order clouds, layer 1).
- R = L L^dag = A A^dag, A = P_Co E H0^{-1/2} on coexact(l'-2q+1) (gapped). So rank R = dim coexact(l'-2q+1) and there
  are no resonant denominators.
Models for R, each inserted into the exact pencil with the exact clouds K:
  exact    R itself
  haar     U R U^dag, U Haar on the edge multiplet (generic position of range R)
  flat     Haar-rotated, rank r_R, all nonzero eigenvalues equal to their mean
  cR_Ew    Haar-rotated diag(c_R / E_w), with c_R the mean of <w|E^dag P_Co E|w> and E_w the exact energies of
           coexact(l'-2q+1)
  diag     diag(R) in the edge basis (D6 per-mode mean field)
  zero     R = 0 (first order)
Full chord model ("chord"), with no finite-N input beyond (n, d_h, r_R):
  - clouds K = S^{1/2} Z / sqrt(d), S from the layer-1 chord population at (s_out, s_in) = (-1/6, 1/2);
  - R = U diag(rho) U^dag, rank r_R (exact count). rho_m = c_R^chord / (kappa_w Lambda_{7/6}(phi_m)), using the
    mirrored coexact family at s = -7/6 at B-family quantiles. kappa_w is calibrated as in D7. The eta-averaged
    c_R^chord = 2 C(N-k, 2) (n / C(N, k+2)), k = degree of coexact(l'-2q+1).
Per (N', seed): T2, H, the fractions with lambda < 0.003, 0.01, 0.03, 0.1, lambda_min, and histograms. Spectra go to .npz.

Usage (under scripts/memwatch.py):
  .venv/bin/python scripts/run_chord_layer3.py --Np 10 12 14 --seeds 9000 9001 9002 --haar 3
  .venv/bin/python scripts/run_chord_layer3.py --Np 16 --seeds 9000 9001 --haar 2
"""
import argparse
import json
import os
import subprocess
import sys
import time
from math import comb, exp

import numpy as np
import scipy.linalg as sla

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from chord_chaos import WedgeStructure, random_couplings  # noqa: E402
from chord_dressing import chord_population  # noqa: E402
from dssyk_n2 import b_family_density, chord_energy  # noqa: E402
from uplift import OldModel  # noqa: E402

OUT = os.path.join(ROOT, "results/data/chord_layer3_2026-10-26")


def git_state():
    try:
        g = subprocess.run(["git", "-C", ROOT, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", ROOT, "status", "--porcelain", "--", "src", "scripts"],
                               capture_output=True, text=True).stdout.strip()
        return g + ("-dirty" if dirty else "")
    except OSError:
        return "unknown"


def haar(n, rng):
    """Haar unitary (QR of a complex Ginibre matrix with phase fix)."""
    Z = (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    d = np.diag(R)
    return Q * (d / np.abs(d))[None, :]


def hblock_lambda(Kv, Rv):
    """Self-dressed h-block decoder eigenvalues from clouds Kv (n x d, edge basis) and R (n x n)."""
    n = Rv.shape[0]
    Iden = np.eye(n)
    A1 = np.linalg.solve(Iden + Rv, Kv)
    G = np.eye(Kv.shape[1]) + Kv.conj().T @ A1
    Nm = np.eye(Kv.shape[1]) - A1.conj().T @ (Iden - Rv) @ A1
    mu = sla.eigvalsh((Nm + Nm.conj().T) / 2, (G + G.conj().T) / 2)
    return np.sort((1 + mu) / 2)


def summary(lam):
    y = (2 * lam - 1) ** 2
    out = dict(T2=float(y.mean()), H=float(np.mean(y * (1 - y))), lam_min=float(lam.min()), mean_lam=float(lam.mean()))
    for x in (0.003, 0.01, 0.03, 0.1):
        out[f"P_lt_{x}"] = float(np.mean(lam < x))
    out["hist"] = np.histogram(lam, bins=40, range=(0, 1))[0].tolist()
    return out


def chord_R_spectrum(N, k, r, cR, nphi=200001):
    """r nonzero eigenvalues rho = cR / (kappa Lambda_{s}(phi)) of the chord back-reaction channel.
    The coexact family at degree k (charge s = (k - N/2)/3) is mirrored to the exact family at -s, as in D7 for the g
    family. phi sits at quantiles of the B-family density, and kappa is calibrated by D7's first-moment rule."""
    p = 3
    lam = 2 * p * p / N
    q = exp(-lam)
    s = -((k - N / 2) / p)                                   # mirrored charge (> 1/2: gapped)
    phi = np.linspace(0, np.pi, nphi)[1:-1]
    dens = b_family_density(lam, s, phi)
    cdf = np.concatenate([[0.0], np.cumsum((dens[1:] + dens[:-1]) / 2 * np.diff(phi))])
    cdf /= cdf[-1]
    ph = np.interp((np.arange(r) + 0.5) / r, cdf, phi)
    deg_in = N - k                                           # mirror degree
    kappa = 2 * (comb(deg_in, p) + comb(N - deg_in, p)) / (q ** s + q ** -s)
    E = kappa * chord_energy(lam, s, ph)
    return cR / E, dict(s_mirror=s, kappa=kappa, lam=lam)


def run(Np, seed, nhaar, rng_seed):
    t0 = time.time()
    rng = np.random.default_rng(seed)
    c = random_couplings(rng, comb(Np, 3))
    F = (Np - 2) // 2
    N = Np - 1
    om = OldModel(Np, c, rng, F=F)
    U, s, Vh = np.linalg.svd(om.Q0_Fm1, full_matrices=False)
    keep = s > 1e-10 * s.max()
    U, s, V = U[:, keep], s[keep], Vh[keep].conj().T
    n = len(s)
    Kv = (U.conj().T @ (om.E_F @ om.VF)) / s[:, None]
    _, _, L = om.three_operators()
    Rv = V.conj().T @ (L @ L.conj().T) @ V
    Rv = (Rv + Rv.conj().T) / 2
    rho = np.linalg.eigvalsh(Rv)
    tol = 1e-10 * rho.max()
    rR = int(np.sum(rho > tol))
    rho_nz = rho[rho > tol]
    # gapped channel: coexact(l'-5) at degree k = F - 3, energies E_w, and c_R
    k = F - 3
    U2, s2, Vh2 = np.linalg.svd(om.Q0_Fm3, full_matrices=False)
    k2 = s2 > 1e-10 * s2.max()
    Wc, Ew = Vh2[k2].conj().T, s2[k2] ** 2
    Emap = WedgeStructure(N, k, 2).matrix(om.eta).toarray()      # L^{F-3} -> L^{F-1}
    EW = Emap @ Wc
    PE = V.conj().T @ EW                                           # P_Co E w in the edge basis
    cR = float(np.mean(np.sum(np.abs(PE) ** 2, axis=0)))
    rec = dict(Np=Np, N=N, F=F, seed=seed, n=n, d_h=int(Kv.shape[1]), rank_R=rR, rank_frac=rR / n,
               rank_frac_chord=float(np.exp(-(10 / 9) * 18 / N)), Emin_rel=float((s ** 2).min() / (s ** 2).mean()),
               cR=cR, cR_chord=float(2 * comb(N - k, 2) * n / comb(N, k + 2)),
               rho_mean=float(rho_nz.mean()), rho_q=[float(x) for x in np.quantile(rho_nz, [0.1, 0.5, 0.9])])
    rng2 = np.random.default_rng(rng_seed)
    models = {}
    models["exact"] = hblock_lambda(Kv, Rv)
    models["zero"] = hblock_lambda(Kv, 0 * Rv)
    models["diag"] = hblock_lambda(Kv, np.diag(np.diag(Rv)))
    flat = np.zeros(n)
    flat[-rR:] = rho_nz.mean()
    mod = np.zeros(n)
    mod[-rR:] = np.sort(cR / Ew)[::-1][:rR] if len(Ew) >= rR else np.pad(np.sort(cR / Ew), (rR - len(Ew), 0))
    rch, info = chord_R_spectrum(N, k, rR, rec["cR_chord"])
    chordR = np.zeros(n)
    chordR[-rR:] = np.sort(rch)
    rec["rho_chord_mean"] = float(rch.mean())
    rec["chord_info"] = info
    # chord clouds: layer-1 population S (n values) with Z (n x d): Kv_chord = S^{1/2} Z / sqrt(d)
    gpop, P = chord_population(Np, "h", n=n, F=F)
    for name in ("haar", "flat", "cR_Ew", "chord"):
        models[name] = []
    for t in range(nhaar):
        Q = haar(n, rng2)
        models["haar"].append(hblock_lambda(Kv, Q @ Rv @ Q.conj().T))
        models["flat"].append(hblock_lambda(Kv, (Q * flat[None, :]) @ Q.conj().T))
        models["cR_Ew"].append(hblock_lambda(Kv, (Q * mod[None, :]) @ Q.conj().T))
        d = Kv.shape[1]
        Z = (rng2.standard_normal((n, d)) + 1j * rng2.standard_normal((n, d))) / np.sqrt(2 * d)
        Kc = np.sqrt(np.clip(gpop, 0, None))[:, None] * Z
        Q2 = haar(n, rng2)
        models["chord"].append(hblock_lambda(Kc, (Q2 * chordR[None, :]) @ Q2.conj().T))
    for name in ("haar", "flat", "cR_Ew", "chord"):
        models[name] = np.sort(np.concatenate(models[name]))
    rec["models"] = {k_: summary(v) for k_, v in models.items()}
    rec["seconds"] = round(time.time() - t0, 1)
    return rec, models


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Np", type=int, nargs="+", required=True)
    ap.add_argument("--seeds", type=int, nargs="+", default=[9000, 9001, 9002])
    ap.add_argument("--haar", type=int, default=3, help="Haar rotations (and chord draws) per realization")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    out = a.out or os.path.join(OUT, "layer3.jsonl")
    meta = dict(argv=sys.argv, git=git_state(), started=time.strftime("%Y-%m-%dT%H:%M:%S"), out=out)
    with open(os.path.join(OUT, f"meta_{time.strftime('%Y%m%dT%H%M%S')}_{os.getpid()}.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    for Np in a.Np:
        for seed in a.seeds:
            rec, models = run(Np, seed, a.haar, rng_seed=31 * seed + Np)
            rec["git"] = meta["git"]
            with open(out, "a") as fh:
                fh.write(json.dumps(rec) + "\n")
            np.savez_compressed(os.path.join(OUT, f"lam_N{Np}_s{seed}.npz"), **models)
            m = rec["models"]
            print(f"N'={Np} seed={seed}: rank R {rec['rank_R']}/{rec['n']} ({rec['rank_frac']:.3f}; chord {rec['rank_frac_chord']:.3f}),"
                  f" c_R {rec['cR']:.1f} (chord {rec['cR_chord']:.1f}), rho mean {rec['rho_mean']:.3f} (chord {rec['rho_chord_mean']:.3f}) [{rec['seconds']}s]",
                  flush=True)
            for name in ("exact", "haar", "flat", "cR_Ew", "chord", "diag", "zero"):
                x = m[name]
                print(f"   {name:6s}: T2 {x['T2']:.4f}  H {x['H']:.4f}  P(<.01) {x['P_lt_0.01']:.3f}  P(<.03) {x['P_lt_0.03']:.3f}"
                      f"  P(<.1) {x['P_lt_0.1']:.3f}  min {x['lam_min']:.2e}", flush=True)


if __name__ == "__main__":
    main()
