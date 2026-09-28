#!/usr/bin/env python3
r"""
Reproduce Table 1 (exact-lift loci) of the metric-fortuity paper for the generic
N=2 SYK model, realized as wedging with a generic q=3 form on the exterior algebra.

Model.  H = Lambda^*(C^N), charge p = degree (dim H^p = C(N,p)).  The supercharge
is Q = C ^ (.) for a generic 3-form C = sum_{i<j<k} C_{ijk} e_i^e_j^e_k, so
Q : H^p -> H^{p+3} and Q^2 = C^C^ = 0 (q=3 odd).  Q^dagger is contraction with Cbar
(the conjugate-transpose in the orthonormal e_S basis).  BPS = ker Q ∩ ker Q^dagger.

Uplift N -> N+1 inherits C and draws independent new couplings C_{ij,N+1}.
  * down channel: decoder reads the e_{N+1}-EMPTY component of B_{N+1}^{p}   (targets ℓ-1)
  * up   channel: decoder reads the e_{N+1}-OCCUPIED component of B_{N+1}^{p+1} (targets ℓ+1)

For each sector the row is
  N  p  ℓ  dim B_N^p | [down: room rkD gen. dimA codim] | [up: room rkD gen. dimA codim]
with room = dim of the SUSY-constrained space the comparison lives in (69),
gen. = dim B + rk D - dim(room) the generic-position prediction, dim A = actual
intersection B ∩ im D, codim = dim B - dim A.

Sectors printed: every (N,p) with ℓ = 2p-N <= 0 and dim B_N^p > 0 (matches the
paper; the empty window-edge sectors, e.g. (7,2), are dropped automatically).

    python reproduce_table1.py --nmin 5 --nmax 6
"""
import argparse, itertools
from math import comb
import numpy as np

Q_BODY = 3

def subsets(n, p): return list(itertools.combinations(range(1, n+1), p))

def merge_sign(T, S):                       # sign of e_T ^ e_S -> e_sorted(T∪S)
    return -1 if (sum(1 for t in T for s in S if t > s) & 1) else 1

def wedge_block(form, n, p):                # matrix of  C ^ (.) : Lambda^p -> Lambda^{p+q}
    cols = subsets(n, p)
    rows = subsets(n, p+Q_BODY)
    rmap = {s: i for i, s in enumerate(rows)}
    M = np.zeros((len(rows), len(cols)), complex)
    for cj, S in enumerate(cols):
        Sset = set(S)
        for T, c in form.items():
            if Sset & set(T): continue
            M[rmap[tuple(sorted(T + S))], cj] += c * merge_sign(T, S)
    return M

def nullspace(M, tol=1e-8):
    if M.shape[0] == 0: return np.eye(M.shape[1], dtype=complex)
    U, s, Vh = np.linalg.svd(M)
    return Vh[int((s > tol).sum()):].conj().T

def rank(M, tol=1e-8):
    if min(M.shape) == 0: return 0
    return int((np.linalg.svd(M, compute_uv=False) > tol).sum())

def harmonic(form, n, p):                    # ker(Q|p) ∩ ker(Q^dagger|p) inside Lambda^p
    d = len(subsets(n, p))
    blocks = []
    Qp = wedge_block(form, n, p)
    if Qp.shape[0]: blocks.append(Qp)
    if p - Q_BODY >= 0:
        blocks.append(wedge_block(form, n, p-Q_BODY).conj().T)
    M = np.vstack(blocks) if blocks else np.zeros((0, d), complex)
    return nullspace(M)

def gen_form(n, seed):
    rng = np.random.default_rng(seed)
    return {T: rng.standard_normal() + 1j*rng.standard_normal()
            for T in itertools.combinations(range(1, n+1), Q_BODY)}

def extend_form(formN, N, seed):             # inherit C, add couplings touching mode N+1
    rng = np.random.default_rng(seed + 777)
    f = dict(formN)
    for combo in itertools.combinations(range(1, N+1), Q_BODY-1):
        f[combo + (N+1,)] = rng.standard_normal() + 1j*rng.standard_normal()
    return f

def component(vecs, N, p, occupied):         # pull e_{N+1}-occupied/empty part into Lambda^p(C^N)
    deg = p+1 if occupied else p
    src = subsets(N+1, deg)
    tmap = {s: i for i, s in enumerate(subsets(N, p))}
    out = np.zeros((len(subsets(N, p)), vecs.shape[1]), complex)
    for i, S in enumerate(src):
        if ((N+1) in S) == occupied:
            key = tuple(x for x in S if x != N+1) if occupied else S
            out[tmap[key], :] += vecs[i, :]  # overall (-1)^p sign is irrelevant to the span
    return out

def channel(B, Bt, N, p, occupied, room_dim):
    D = component(Bt, N, p, occupied)                 # im D in Lambda^p(C^N)
    dimB, rkD = B.shape[1], rank(D)
    dimBplus = rank(np.hstack([B, D])) if D.shape[1] else dimB
    dimA = dimB + rkD - dimBplus
    gen  = max(dimB + rkD - room_dim, 0)
    return (room_dim, rkD, gen, dimA, dimB - dimA)

def row(N, p, seed):
    fN  = gen_form(N, seed)
    B   = harmonic(fN, N, p)
    dimB = B.shape[1]
    if dimB == 0: return None
    fN1 = extend_form(fN, N, seed)
    rkMp = rank(wedge_block(fN, N, p))
    down_room = comb(N, p) - rkMp                     # dim ker Q_N|p
    up_room   = dimB + rkMp                            # dim ker Q_N^dagger|p
    Bd = harmonic(fN1, N+1, p)                          # target for down channel
    Bu = harmonic(fN1, N+1, p+1)                        # target for up channel
    return {'N': N, 'p': p, 'ell': 2*p - N, 'dimB': dimB,
            'down': channel(B, Bd, N, p, False, down_room),
            'up':   channel(B, Bu, N, p, True,  up_room)}

def print_table(nmin, nmax, seed):
    C = ['room', 'rk D', 'gen.', 'dim A', 'codim']
    w = 6
    top = (f"{'':>3} {'':>3} {'':>4} {'':>9} | "
           f"{'down channel':^{5*(w+1)-1}} | {'up channel':^{5*(w+1)-1}}")
    hdr = (f"{'N':>3} {'p':>3} {'ℓ':>4} {'dim B_N^p':>9} | "
           + " ".join(f"{c:>{w}}" for c in C) + " | "
           + " ".join(f"{c:>{w}}" for c in C))
    print(top); print(hdr); print('-'*len(hdr))
    for N in range(nmin, nmax+1):
        for p in range(0, N//2 + 1):                   # ℓ = 2p-N <= 0
            r = row(N, p, seed)
            if r is None: continue                     # drop empty BPS sectors
            d, u = r['down'], r['up']
            print(f"{r['N']:>3} {r['p']:>3} {r['ell']:>4} {r['dimB']:>9} | "
                  + " ".join(f"{x:>{w}}" for x in d) + " | "
                  + " ".join(f"{x:>{w}}" for x in u))

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nmin", type=int, default=5)
    ap.add_argument("--nmax", type=int, default=6)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    print_table(a.nmin, a.nmax, a.seed)
