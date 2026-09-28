#!/usr/bin/env python3
r"""
protected_atoms.py -- exact count of the decoder's protected atoms (lambda=0,1)
for the one-flavor N=2 SYK uplift decoder N -> N+1, as a function of N.

WHAT IT COMPUTES.  The exactly-transmitted (lambda=1) and exactly-lost (lambda=0)
BPS directions of the uplift decoder are, in closed operator form,

    #{lambda=1} = dim( B_N^p  cap  ker(D wedge : Lam^p    -> Lam^{p+2}) ),
    #{lambda=0} = dim( B_N^{p-1} cap ker(iota_D : Lam^{p-1} -> Lam^{p-3}) ),

with p = (N+1)//2, where
  * B_N^k = harmonic k-forms of the OLD (size-N) theory: ker(C wedge) cap ker(C wedge)^dagger,
    C a generic complex 3-form (the inherited couplings);
  * D is a generic complex 2-form (the new couplings C_{ij,N+1}); iota_D = (D wedge)^dagger.
These reproduce the exact-diagonalization atom counts (verified: N=12 -> 16,25 ; N=13 -> 41,41).

WHY.  We want the large-N fate of the protected fraction (atoms)/d, where the enlarged
BPS dimension is d = 3^{N/2} (N even), 2*3^{(N-1)/2} (N odd).  Finding so far:
  N=12: 5.62%   N=13: 5.62%   N=14: 0.00%  (two seeds)
i.e. the atoms look like a finite-N resonance (peaking where dim B_N^p ~ C(N,p+2), around
N=12-13) that vanishes for N>=14.  THIS SCRIPT IS TO CONFIRM N=15,16.

METHOD.  Build the harmonic space B_N^k as the nullspace of the sparse constraint
[ C^:k->k+3 ; (C^:k-3->k)^dagger ], then the atom count is the nullity of the new-coupling
map restricted to it, obtained from the eigenvalues of the small (dimB x dimB) Gram matrix.

COST (dense backend).  Dominant memory is the nullspace SVD, whose V factor is
C(N,p) x C(N,p) complex (16 bytes/entry):
    N=15: C(15,8)=6435   -> ~0.66 GB,  ~1-2 min
    N=16: C(16,8)=12870  -> ~2.65 GB,  ~5-10 min
For N=16 on a memory-tight machine use --backend sparseqr (pip install sparseqr), which
avoids forming the full V; verify it against dense at N=12 first.

USAGE
    python protected_atoms.py --N 15 16 --seeds 2
    python protected_atoms.py --N 16 --seeds 1 --backend sparseqr
    python protected_atoms.py --N 12 13 --seeds 1          # reproduces 16,25 / 41,41
"""
import argparse, itertools, os, time
from math import comb
import numpy as np
import scipy.sparse as sp

Q = 3  # supercharge is a 3-form


# ----------------------------------------------------------------- combinatorics
def subset_index(N, p):
    masks = []
    for S in itertools.combinations(range(N), p):
        m = 0
        for i in S:
            m |= (1 << i)
        masks.append(m)
    return masks, {m: i for i, m in enumerate(masks)}


def merge_sign(Tmask, Smask):
    """Sign of e_T ^ e_S -> e_{sorted(T u S)} for disjoint T, S (bitmasks)."""
    s = 0
    t = Tmask
    while t:
        b = t & (-t)
        s += (Smask & (b - 1)).bit_count()
        t ^= b
    return -1.0 if (s & 1) else 1.0


def form_wedge(coeffs, N, p, deg):
    """Sparse matrix of (omega ^ .) : Lam^p -> Lam^{p+deg}, omega = {tuple(modes): value}."""
    if p < 0 or p + deg > N:
        return sp.csr_matrix((comb(N, p + deg) if 0 <= p + deg <= N else 0,
                              comb(N, p) if 0 <= p <= N else 0), dtype=complex)
    cols, _ = subset_index(N, p)
    _, rmap = subset_index(N, p + deg)
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
                I.append(rmap[Smask | Tmask])
                J.append(cj)
                V.append(c * merge_sign(Tmask, Smask))
    return sp.csr_matrix((V, (I, J)), shape=(comb(N, p + deg), comb(N, p)), dtype=complex)


def rand_form(seed, N, deg):
    rng = np.random.default_rng(seed)
    d = {}
    for T in itertools.combinations(range(N), deg):
        d[T] = rng.standard_normal() + 1j * rng.standard_normal()
    return d


# ----------------------------------------------------------------- linear algebra
def nullspace_dense(M):
    A = M.toarray()
    if min(A.shape) == 0:
        return np.eye(A.shape[1], dtype=complex)
    U, s, Vh = np.linalg.svd(A, full_matrices=True)
    tol = max(A.shape) * (s[0] if len(s) else 1.0) * np.finfo(float).eps
    rank = int((s > tol).sum())
    return Vh[rank:].conj().T


def nullspace_sparseqr(M, tol=1e-9):
    import sparseqr
    Mt = M.conj().transpose().tocoo()
    Qm, R, E, rank = sparseqr.qr(Mt)
    return Qm.toarray()[:, rank:]


def harmonic(C, N, p, backend="dense"):
    """Orthonormal basis (columns) of B_N^p = ker(C^) cap ker(C^)^dagger at degree p."""
    if p < 0 or p > N:
        return np.zeros((max(comb(N, p) if 0 <= p <= N else 0, 0), 0), dtype=complex)
    blocks = [form_wedge(C, N, p, Q)]
    if p - Q >= 0:
        blocks.append(form_wedge(C, N, p - Q, Q).conj().transpose())
    M = sp.vstack(blocks).tocsr()
    B = (nullspace_sparseqr if backend == "sparseqr" else nullspace_dense)(M)
    if B.shape[1]:
        B, _ = np.linalg.qr(B)
    return B


def nullity_on(basis, op):
    """dim( span(basis) cap ker(op) ) via eigenvalues of the (dimB x dimB) Gram matrix."""
    if basis.shape[1] == 0:
        return 0
    X = op @ basis                          # (target x dimB)
    G = X.conj().T @ X                       # (dimB x dimB), Hermitian PSD
    ev = np.linalg.eigvalsh((G + G.conj().T) / 2).real
    tol = max(G.shape) * (ev[-1] if len(ev) else 1.0) * np.finfo(float).eps * 1e3
    rank = int((ev > tol).sum())
    return basis.shape[1] - rank


# ----------------------------------------------------------------- counts
def d_enl(N):
    return 3 ** (N // 2) * (2 if N % 2 else 1)


def bdim(N, k):
    def Dtil(M, kk):
        s, n = 0, 0
        while kk - n * Q >= 0:
            s += (-1) ** n * comb(M, kk - n * Q)
            n += 1
        return s
    if k < 0 or k > N:
        return 0
    return Dtil(N, k) - Dtil(N, N - k - Q)


def atom_counts(N, seed, backend="dense"):
    p = (N + 1) // 2
    C = rand_form(seed, N, 3)
    D = rand_form(seed + 100003, N, 2)             # independent new couplings
    Bp = harmonic(C, N, p, backend)                # B_N^p
    lam1 = nullity_on(Bp, form_wedge(D, N, p, 2))  # ker(D^: p -> p+2)
    Bpm = harmonic(C, N, p - 1, backend)           # B_N^{p-1}
    iota = form_wedge(D, N, p - 3, 2).conj().transpose()  # iota_D: p-1 -> p-3
    lam0 = nullity_on(Bpm, iota)
    return lam1, lam0, Bp.shape[1], Bpm.shape[1]


# ----------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--N", type=int, nargs="+", required=True)
    ap.add_argument("--seeds", type=int, default=1)
    ap.add_argument("--backend", choices=["dense", "sparseqr"], default="dense")
    ap.add_argument("--out", default="protected_atoms.csv")
    a = ap.parse_args()
    if not os.path.exists(a.out):
        open(a.out, "w").write("N,seed,p,lam1,lam0,atoms,d,frac,dimBp,dimBpm1,secs\n")
    print(f"{'N':>3} {'seed':>4} {'p':>3} {'lam=1':>6} {'lam=0':>6} {'atoms':>6} "
          f"{'d':>8} {'frac':>7} | {'dimB^p':>7} {'dimB^p-1':>8} {'C(N,p+2)':>9}")
    for N in a.N:
        p = (N + 1) // 2
        for s in range(a.seeds):
            t0 = time.time()
            l1, l0, bp, bpm = atom_counts(N, s, a.backend)
            secs = time.time() - t0
            d = d_enl(N); tot = l1 + l0
            print(f"{N:>3} {s:>4} {p:>3} {l1:>6} {l0:>6} {tot:>6} {d:>8} "
                  f"{tot/d:>6.2%} | {bp:>7} {bpm:>8} {comb(N, p+2) if p+2<=N else 0:>9}"
                  f"   ({secs:.0f}s)")
            with open(a.out, "a") as f:
                f.write(f"{N},{s},{p},{l1},{l0},{tot},{d},{tot/d:.5f},{bp},{bpm},{secs:.1f}\n")


if __name__ == "__main__":
    main()
