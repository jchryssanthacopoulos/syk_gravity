#!/usr/bin/env python3
r"""
decoder_fast.py -- optimized builder for the one-flavor N=2 SYK uplift decoder,
to push the soft-tail study to larger N.

MODEL.  Enlarged theory N' = N+1, exterior algebra H = Lambda^*(C^{N'}), supercharge
Q = C ^ (.) with a generic 3-form C.  BPS space B = ker Q cap ker Q^dagger at degree
P = (N+1)//2 (half filling of the enlarged theory).  The down-channel decoder is

        O = P_B (1 - n_{N'}) P_B     restricted to B,        eigenvalues lambda in [0,1],

where (1 - n_{N'}) projects onto "new mode N' empty".  lambda = <psi|(1-n_{N'})|psi>.

WHY THIS IS FASTER THAN THE NAIVE DENSE BUILD.
  * The wedge operator Q is built sparsely with integer bitmasks + popcount signs
    (the pure-Python dense double loop was the real bottleneck).
  * The BPS space is the nullspace of the SPARSE constraint matrix
        M = [ Q_P ;  Q_{P-3}^dagger ]              (ker Q cap ker Q^dagger),
    never the C(N',P)^2 Hodge Laplacian.
  * The decoder eigenvalues are the SQUARED SINGULAR VALUES of the "N'-empty" block
    B_empty of the orthonormal BPS basis B.  Using singular values (not the Gram
    matrix B_empty^dagger B_empty) keeps full float64 precision on the tiny tail
    eigenvalues (sigma ~ sqrt(lambda) ~ 1e-3 stays well away from 1e-16).

NULLSPACE BACKENDS (--backend).
  'dense'    scipy dense SVD of the sparse-built M.  Portable.  Memory ~ 16*C(N',P)^2
             bytes (the SVD's V factor): ~0.7 GB at N'=15, ~2.6 GB at N'=16,
             ~9 GB at N'=17.  This is the recommended default up to N'=16.
  'sparseqr' SuiteSparseQR (pip install sparseqr).  Much lower memory; use for
             N' >= 17.  See nullspace_sparseqr() -- verify against 'dense' at a small
             N' first, since the exact idiom can vary by package version.

USAGE
  python decoder_fast.py --N 12 13 14 --seeds 3
  python decoder_fast.py --N 15 --seeds 1 --backend dense --out spec
Outputs spec_N{N}_s{seed}.npy (sorted eigenvalues) and appends a line per run to
softtail_scaling.csv with the tail statistics.
"""
import argparse, itertools, time, os
from math import comb, sqrt, pi
import numpy as np
import scipy.sparse as sp

q = 3  # supercharge is a 3-form

# ------------------------------------------------------------------ combinatorics
def subset_index(Nprime, p):
    """Return (list of p-subset bitmasks, dict bitmask->row index) in colex order."""
    masks = []
    for S in itertools.combinations(range(Nprime), p):
        m = 0
        for i in S:
            m |= (1 << i)
        masks.append(m)
    return masks, {m: i for i, m in enumerate(masks)}

def merge_sign(Tmask, Smask):
    """Sign of e_T ^ e_S -> e_{sorted(T u S)} for disjoint T,S (bitmasks)."""
    s = 0
    t = Tmask
    while t:
        b = t & (-t)              # lowest set bit of T
        i = b.bit_length() - 1
        s += (Smask & (b - 1)).bit_count()   # # of S-elements below position i
        t ^= b
    return -1.0 if (s & 1) else 1.0

def couplings(seed, Nprime):
    """Generic complex 3-form: dict {3-subset bitmask: C_ijk}."""
    rng = np.random.default_rng(seed)
    C = {}
    for T in itertools.combinations(range(Nprime), q):
        m = 0
        for i in T:
            m |= (1 << i)
        C[m] = rng.standard_normal() + 1j * rng.standard_normal()
    return C

# ------------------------------------------------------------------ wedge operator
def wedge(C, Nprime, p):
    """Sparse matrix of  C ^ (.) : Lambda^p -> Lambda^{p+q},  shape C(N',p+q) x C(N',p)."""
    cols, _ = subset_index(Nprime, p)
    _, rmap = subset_index(Nprime, p + q)
    I, J, V = [], [], []
    for cj, S in enumerate(cols):
        for T, c in C.items():
            if S & T:                      # T must be disjoint from S
                continue
            I.append(rmap[S | T]); J.append(cj); V.append(c * merge_sign(T, S))
    shape = (comb(Nprime, p + q), comb(Nprime, p))
    return sp.csr_matrix((V, (I, J)), shape=shape, dtype=complex)

# ------------------------------------------------------------------ BPS nullspace
def nullspace_dense(M, tol=None):
    """Orthonormal nullspace of sparse M via dense SVD (columns)."""
    A = M.toarray()
    # full SVD; nullspace = right-singular vectors with ~zero singular value
    U, s, Vh = np.linalg.svd(A, full_matrices=True)
    if tol is None:
        tol = max(A.shape) * (s[0] if len(s) else 1.0) * np.finfo(float).eps
    rank = int((s > tol).sum())
    return Vh[rank:].conj().T          # (cols x nullity)

def nullspace_sparseqr(M, tol=1e-9):
    """Orthonormal nullspace via SuiteSparseQR (advanced; verify vs dense once)."""
    import sparseqr
    # QR of M^H = Q R ;  columns of Q past the rank span ker(M).
    Mt = M.conj().transpose().tocoo()
    Q, R, E, rank = sparseqr.qr(Mt)
    Qd = Q.toarray()
    return Qd[:, rank:]

def bps_basis(C, Nprime, P, backend="dense"):
    """Orthonormal basis (cols) of ker Q_P cap ker Q^dagger_P at degree P."""
    blocks = [wedge(C, Nprime, P)]                       # Q_P : Lambda^P -> Lambda^{P+q}
    if P - q >= 0:
        blocks.append(wedge(C, Nprime, P - q).conj().transpose())  # Q^dagger_P
    M = sp.vstack(blocks).tocsr()
    B = (nullspace_sparseqr if backend == "sparseqr" else nullspace_dense)(M)
    # re-orthonormalize (cheap; guards against backend tolerance)
    Bq, _ = np.linalg.qr(B)
    return Bq

# ------------------------------------------------------------------ decoder + stats
def Dtil(Mval, k):
    s, n = 0, 0
    while k - n * q >= 0:
        s += (-1) ** n * comb(Mval, k - n * q); n += 1
    return s

def decoder_spectrum(Nprime, P, seed, backend="dense"):
    """Sorted decoder eigenvalues lambda = sigma(B_empty)^2."""
    C = couplings(seed, Nprime)
    B = bps_basis(C, Nprime, P, backend)                 # (C(N',P) x d)
    d = B.shape[1]
    masks, _ = subset_index(Nprime, P)
    top = 1 << (Nprime - 1)
    empty_rows = [i for i, m in enumerate(masks) if not (m & top)]   # mode N'-1 absent
    Be = B[empty_rows, :]                                # (n_empty x d)
    sig = np.linalg.svd(Be, compute_uv=False)            # length min(n_empty, d)
    lam = np.sort(np.clip(sig, 0, 1) ** 2)
    if len(lam) < d:                                     # occupied-only -> exact lambda=0
        lam = np.concatenate([np.zeros(d - len(lam)), lam])
    return lam, d

def wachter_cdf(a, b, xs):
    """Normalized CDF of the interior Wachter density on [lam_-, lam_+]."""
    lm = (sqrt(a * (1 - b)) - sqrt(b * (1 - a))) ** 2
    lp = (sqrt(a * (1 - b)) + sqrt(b * (1 - a))) ** 2
    g = np.linspace(lm, lp, 6000)[1:-1]
    f = np.sqrt(np.clip((lp - g) * (g - lm), 0, None)) / (2 * pi * a * g * (1 - g))
    cdf = np.concatenate([[0.0], np.cumsum((f[1:] + f[:-1]) / 2 * np.diff(g))])
    cdf /= cdf[-1]
    return np.interp(xs, g, cdf, left=0.0, right=1.0)

def tail_stats(lam, Nprime, P, d, tol0=1e-9):
    N = Nprime - 1
    D = comb(Nprime, P); m = comb(N, P)
    a = d / D; b = m / D
    lm = (sqrt(a * (1 - b)) - sqrt(b * (1 - a))) ** 2
    lp = (sqrt(a * (1 - b)) + sqrt(b * (1 - a))) ** 2
    gen = lam[lam >= tol0]
    tail = gen[gen < lm]
    kappa = float(np.mean(1.0 / gen) - 1.0) if len(gen) else float("nan")
    # (1) deep-tail fraction: below lam_-/5  (lambda_- -independent of small shifts)
    frac_deep = float(np.sum(gen < lm / 5) / d)
    # (2) KS of the BULK (interior [lam_-,lam_+]) to Wachter -- should stay small if the
    #     bulk is universal while the tail carries the deviation
    interior = np.sort(gen[(gen >= lm) & (gen <= lp)])
    ks = float("nan")
    if len(interior) > 5:
        Fe = np.arange(1, len(interior) + 1) / len(interior)
        ks = float(np.max(np.abs(Fe - wachter_cdf(a, b, interior))))
    # tail density exponent: cumulative count C(l) ~ l^beta => rho ~ l^{beta-1}
    beta = float("nan")
    if len(tail) > 12:
        ts = np.sort(tail); Cc = np.arange(1, len(ts) + 1)
        lo, hi = 3, len(ts) - 2
        beta = float(np.polyfit(np.log(ts[lo:hi]), np.log(Cc[lo:hi]), 1)[0])
    return dict(N=N, d=d, a=a, b=b, lam_min=(gen[0] if len(gen) else np.nan),
                lam_edge=lm, n_tail=int(len(tail)), tail_frac=len(tail) / d,
                frac_deep=frac_deep, ks_bulk=ks,
                kappa=kappa, rho_exp=beta - 1 if beta == beta else np.nan)

# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--N", type=int, nargs="+", required=True, help="values of N (uplift N->N+1)")
    ap.add_argument("--seeds", type=int, default=1, help="realizations per N")
    ap.add_argument("--backend", choices=["dense", "sparseqr"], default="dense")
    ap.add_argument("--out", default="spec", help="prefix for saved spectra")
    a = ap.parse_args()
    csv = "softtail_scaling.csv"
    if not os.path.exists(csv):
        open(csv, "w").write("N,seed,d,a,b,lam_min,lam_edge,n_tail,tail_frac,"
                             "frac_deep,ks_bulk,kappa,rho_exp,secs\n")
    for N in a.N:
        Nprime, P = N + 1, (N + 1) // 2
        for s in range(a.seeds):
            t0 = time.time()
            lam, d = decoder_spectrum(Nprime, P, s, a.backend)
            np.save(f"{a.out}_N{N}_s{s}.npy", lam)
            st = tail_stats(lam, Nprime, P, d)
            secs = time.time() - t0
            print(f"N={N} seed={s}: d={d} a={st['a']:.3f} b={st['b']:.3f} "
                  f"lam_edge={st['lam_edge']:.4f} lam_min={st['lam_min']:.2e} "
                  f"tail_frac={st['tail_frac']:.4f} frac_deep={st['frac_deep']:.4f} "
                  f"ks_bulk={st['ks_bulk']:.3f} kappa={st['kappa']:.1f} "
                  f"rho_exp={st['rho_exp']:.3f} ({secs:.0f}s)")
            with open(csv, "a") as f:
                f.write(f"{N},{s},{d},{st['a']:.5f},{st['b']:.5f},{st['lam_min']:.6e},"
                        f"{st['lam_edge']:.6f},{st['n_tail']},{st['tail_frac']:.5f},"
                        f"{st['frac_deep']:.5f},{st['ks_bulk']:.5f},"
                        f"{st['kappa']:.4f},{st['rho_exp']:.4f},{secs:.1f}\n")

if __name__ == "__main__":
    main()
