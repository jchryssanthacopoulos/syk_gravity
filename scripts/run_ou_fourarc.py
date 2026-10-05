#!/usr/bin/env python3
r"""
Phase A1/A2 (docs/derivations.md D5): principal-angle spectra of OU-correlated BPS projector pairs, and of the
decoder / parity-family pairs (P, U_A P U_A) built on the *same* realizations.

OU family.  Angle trick (as in run_chord_chaos.py): per sample draw omega, omega'; theta_j = j pi / K; P_j = BPS
projector of omega_j = cos(theta_j) omega + sin(theta_j) omega'.  Every unordered pair (i, j) is an OU pair with
x = |cos((j - i) pi / K)| (P(-omega) = P(omega)).  For each pair the cosines^2 of the principal angles are the
squared singular values of O = V_i^dag V_j (V = orthonormal BPS basis, D x d), so
        T_2m(x) = E tau((P_i P_j)^m) = E mean(c^(2m)),      tau = Tr / d.
Lag K/2 gives x = 0 exactly (independent pair: the two-model value c0); lag 0 gives T = 1.

Site-local family (--family site, D5 Sec. 1 "T(1, z)").  Only the couplings J containing site --site rotate:
C_B(theta) = cos(theta) C_B + sin(theta) C'_B, C_A fixed, theta_j = 2 pi j / K.  A pair at lag m has
z = cos(2 pi m / K) in [-1, 1] (no z -> -z symmetry; z = -1 at m = K/2 IS the decoder of that site, C_B -> -C_B).
Classes are lags m = 1..K/2 (m and K - m give the same z).

Decoder.  On the first --probe-proj projectors of each sample, for U_A = prod_{i in A} (-1)^{n_i}:
O = V^dag U_A V is Hermitian with eigenvalues mu (the spectrum of P U_A P on the BPS space), and
T_2m = mean(mu^(2m)).  A = every single site (the decoder proper), plus --parity-sets random sets for each
|A| = 2 .. N'//2 (the parity family; x_exact(s) = Krawtchouk_p(s) / binom(N', s)).

Recorded per sample (JSON line): a, d; per |x| class: x, number of pairs, mean moments m_k = mean(c^(2k)), k = 1..kmax,
and a histogram of c^2 (summed over pairs); per probe: |A|, sites, moments of mu^2, histogram of mu^2.
Analysis: scripts/analyze_ou_fourarc.py.

Memory: K thin bases (K * D * d * 16 B) + one dense projector and Hodge SVD workspace (~6 D^2 * 16 B).
  N'=12: ~0.1 GB;  N'=14 (K=16): ~2.4 GB.  Use --dry-run, and run under scripts/memwatch.py:
  .venv/bin/python scripts/memwatch.py --limit-gb 4 --log results/data/ou_fourarc_2026-10-05/memwatch.jsonl -- \
      .venv/bin/python scripts/run_ou_fourarc.py --Np 14 --K 16 --samples 6 --seed0 0
"""
import argparse
import json
import os
import subprocess
import sys
import time
from math import comb

import numpy as np
import scipy.linalg as sla

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from chord_chaos import Sector, bps_projector_hodge, p_subsets, random_couplings  # noqa: E402


def git_state():
    try:
        g = subprocess.run(["git", "-C", ROOT, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", ROOT, "status", "--porcelain"], capture_output=True, text=True).stdout.strip()
        return g + ("-dirty" if dirty else "")
    except OSError:
        return "unknown"


def bps_dim(N, p, F):
    """Exact BPS dimension from the cohomology count (BBNR eq. 5.25): sum_n (-1)^n binom(N, F + n p)."""
    return sum((-1) ** abs(n) * comb(N, F + n * p) for n in range(-(F // p) - 1, (N - F) // p + 2)
               if 0 <= F + n * p <= N)


def mem_estimate_gb(N, K, p=3, F=None):
    F = N // 2 if F is None else F
    D = comb(N, F)
    d = bps_dim(N, p, F)
    return (K * D * d * 16 + 6 * D * D * 16 + 2 * D * d * 16) / 1e9 + 0.3


def bps_basis(S, c, d_expected, rng):
    """Orthonormal basis V (D x d) of ker H via the Hodge projector and a randomized range finder (P has exact rank d,
    so range(P G) = range(P)); one re-projection pass.  Returns V and the max deviation ||P V - V||."""
    P, _, d = bps_projector_hodge(S, c, d_expected)
    G = rng.standard_normal((S.D, d)) + 1j * rng.standard_normal((S.D, d))
    V, _ = np.linalg.qr(P @ G)
    V, _ = np.linalg.qr(P @ V)
    err = np.abs(P @ V - V).max()
    del P, G
    return V, err, d


def _orth_range_gram(A, side, rel_tol=1e-10):
    """Orthonormal basis of the column space of A (side='col', D x r) or of A^dag (side='row', A is r x D) via the
    eigendecomposition of the small Gram matrix (BLAS3), followed by a QR re-orthonormalization."""
    G = A @ A.conj().T if side == "row" else A.conj().T @ A
    w, U = sla.eigh(G, check_finite=False)
    keep = w > rel_tol * max(w.max(), 1.0)
    W = U[:, keep] / np.sqrt(w[keep])
    R = A.conj().T @ W if side == "row" else A @ W
    R, _ = np.linalg.qr(R)
    return R, w


def bps_basis_fast(S, c, d_expected, rng):
    """BPS basis without the dense D x D projector or thin SVDs: im Q_F^dag and im Q_{F-p} from Gram eigh (BLAS3),
    then V = orthonormal complement via a randomized projection G - R (R^dag G), two passes.  Returns V, max|Q V|."""
    Qu = S.up.matrix(c).toarray()                        # C(N,F+p) x D
    R1, _ = _orth_range_gram(Qu, "row")
    Rs = [R1]
    Qd = None
    if S.down.shape[0]:
        Qd = S.down.matrix(c).toarray()                  # D x C(N,F-p)
        R2, _ = _orth_range_gram(Qd, "col")
        Rs.append(R2)
    R = np.hstack(Rs)
    d = S.D - R.shape[1]
    if d_expected is not None and d != d_expected:
        raise RuntimeError(f"BPS dimension {d} != expected {d_expected}")
    V = rng.standard_normal((S.D, d)) + 1j * rng.standard_normal((S.D, d))
    for _ in range(2):
        V = V - R @ (R.conj().T @ V)
        V, _ = np.linalg.qr(V)
    err = np.abs(Qu @ V).max()
    if Qd is not None:
        err = max(err, np.abs(Qd.conj().T @ V).max())
    return V, float(err / max(np.abs(Qu).max(), 1.0)), d


def moments_hist(c2, kmax, bins):
    c2 = np.clip(c2, 0.0, 1.0)
    mom = [float(np.mean(c2 ** k)) for k in range(1, kmax + 1)]
    h, _ = np.histogram(c2, bins=bins, range=(0.0, 1.0))
    return mom, h


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Np", type=int, required=True, help="number of complex modes N'")
    ap.add_argument("--p", type=int, default=3)
    ap.add_argument("--F", type=int, default=None, help="charge sector (default N'//2)")
    ap.add_argument("--K", type=int, default=16, help="angles per sample (even; lag K/2 is x = 0)")
    ap.add_argument("--samples", type=int, default=4)
    ap.add_argument("--seed0", type=int, default=0, help="sample s uses numpy seed 1_000_003 * (seed0 + s) + 29")
    ap.add_argument("--probe-proj", type=int, default=2, help="projectors per sample probed with parities")
    ap.add_argument("--parity-sets", type=int, default=3, help="random sets per |A| = 2..N'//2")
    ap.add_argument("--family", choices=["ou", "site"], default="ou")
    ap.add_argument("--site", type=int, default=None, help="site for --family site (default N'-1)")
    ap.add_argument("--fast-basis", action="store_true", help="Gram-eigh basis (BLAS3; for N' >= 14)")
    ap.add_argument("--pair-starts", type=int, default=0,
                    help="if > 0, use only pairs (i, i+m) with i in K/pair_starts * range(pair_starts) (fewer SVDs)")
    ap.add_argument("--kmax", type=int, default=8)
    ap.add_argument("--bins", type=int, default=200)
    ap.add_argument("--out", default=os.path.join(ROOT, "results/data/ou_fourarc_2026-10-05"))
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    F = a.Np // 2 if a.F is None else a.F
    est = mem_estimate_gb(a.Np, a.K, a.p, F)
    print(f"N'={a.Np} p={a.p} F={F} K={a.K}: estimated peak memory {est:.2f} GB", flush=True)
    if a.dry_run:
        return
    assert a.K % 2 == 0
    os.makedirs(a.out, exist_ok=True)
    S = Sector(a.Np, a.p, F=F)
    d_exp = bps_dim(a.Np, a.p, F)
    started = time.strftime("%Y%m%dT%H%M%S")
    tag = f"p{a.p}_N{a.Np}_F{F}" + ("" if a.family == "ou" else f"_site{a.Np - 1 if a.site is None else a.site}")
    meta = dict(args=vars(a), git=git_state(), started=started, D=S.D, d=d_exp, M=S.M, F=F, mem_estimate_gb=est,
                description="OU-pair and decoder principal-angle spectra on matched realizations; docs/derivations.md D5")
    with open(os.path.join(a.out, f"meta_{tag}_{started}.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    masks = np.array(S.masks, dtype=np.int64)
    K = a.K
    if a.family == "ou":
        thetas = np.pi * np.arange(K) / K
    else:
        site = a.Np - 1 if a.site is None else a.site
        _, tm = p_subsets(a.Np, a.p)
        Bmask = np.array([(t >> site) & 1 for t in tm], dtype=bool)
        thetas = 2 * np.pi * np.arange(K) / K
    out_path = os.path.join(a.out, f"ou_{tag}.jsonl")
    for s in range(a.samples):
        t0 = time.time()
        seed = 1_000_003 * (a.seed0 + s) + 29
        rng = np.random.default_rng(seed)
        c1, c2 = random_couplings(rng, S.M), random_couplings(rng, S.M)
        Vs, errs = [], []
        for th in thetas:
            if a.family == "ou":
                c = np.cos(th) * c1 + np.sin(th) * c2
            else:
                c = c1.copy()
                c[Bmask] = np.cos(th) * c1[Bmask] + np.sin(th) * c2[Bmask]
            V, err, d = (bps_basis_fast if a.fast_basis else bps_basis)(S, c, d_exp, rng)
            Vs.append(V)
            errs.append(float(err))
        # OU pairs, grouped by lag class L = min(m, K - m), x = cos(L pi / K)
        cls = {L: dict(n=0, mom=np.zeros(a.kmax), hist=np.zeros(a.bins, dtype=np.int64)) for L in range(1, K // 2 + 1)}
        if a.pair_starts > 0:
            starts = [k * (K // a.pair_starts) for k in range(a.pair_starts)]
            pairs = sorted({(min(i, (i + m) % K), max(i, (i + m) % K)) for i in starts for m in range(1, K // 2 + 1)})
        else:
            pairs = [(i, j) for i in range(K) for j in range(i + 1, K)]
        for i, j in pairs:
            if True:
                m = j - i
                L = min(m, K - m)
                O = Vs[i].conj().T @ Vs[j]
                c2 = sla.eigvalsh(O @ O.conj().T, check_finite=False)   # = squared singular values (faster)
                mom, h = moments_hist(c2, a.kmax, a.bins)
                cls[L]["n"] += 1
                cls[L]["mom"] += np.array(mom)
                cls[L]["hist"] += h
        ang = np.pi / K if a.family == "ou" else 2 * np.pi / K
        ou = [dict(lag=L, x=float(np.cos(L * ang)), pairs=c["n"], moments=(c["mom"] / c["n"]).tolist(),
                   hist=c["hist"].tolist()) for L, c in sorted(cls.items())]
        # decoder / parity probes on the first probe_proj projectors
        probes = []
        prng = np.random.default_rng(seed + 7)
        for jp in range(min(a.probe_proj, K)):
            V = Vs[jp * (K // max(a.probe_proj, 1))]          # spread over the circle (independent-ish)
            sets = [[i] for i in range(a.Np)]
            for size in range(2, a.Np // 2 + 1):
                for _ in range(a.parity_sets):
                    sets.append(sorted(prng.choice(a.Np, size=size, replace=False).tolist()))
            for A in sets:
                Amask = 0
                for i in A:
                    Amask |= (1 << i)
                sign = 1.0 - 2.0 * (np.array([bin(int(x)).count("1") for x in (masks & Amask)]) % 2)
                O = V.conj().T @ (sign[:, None] * V)
                mu = sla.eigvalsh(O, check_finite=False)
                mom, h = moments_hist(mu ** 2, a.kmax, a.bins)
                probes.append(dict(proj=jp, A=A, s=len(A), mu_mean=float(mu.mean()), moments=mom, hist=h.tolist()))
        rec = dict(family=a.family, sample=s, seed=seed, Np=a.Np, p=a.p, F=F, D=S.D, d=d_exp, a=d_exp / S.D, K=K,
                   basis_err=max(errs), ou=ou, probes=probes, secs=time.time() - t0)
        with open(out_path, "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        t2 = {o["lag"]: o["moments"][0] for o in ou}
        t4 = {o["lag"]: o["moments"][1] for o in ou}
        dec = [p_["moments"] for p_ in probes if p_["s"] == 1]
        print(f"sample {s}: {time.time() - t0:.1f}s  basis_err={max(errs):.1e}  T2(x=0)={t2[K // 2]:.4f} "
              f"T4(x=0)={t4[K // 2]:.4f}  decoder T2={np.mean([x[0] for x in dec]):.4f} "
              f"T4={np.mean([x[1] for x in dec]):.4f}", flush=True)
        del Vs


if __name__ == "__main__":
    main()
