r"""
size_decomposition.py -- U(N) "operator-size" decomposition of operators on the charge-P sector Lambda^P(C^N).

Facts used (research/notes/decoder_moments_2026-09-29.md):
* U(N) acts on single-particle modes; on Lambda^P it is irreducible, and End(Lambda^P) = Lambda^P (x) (Lambda^P)^*
  decomposes multiplicity-free into irreps V_k, k = 0..min(P, N-P), with highest weight (1^k, 0, ..., 0, (-1)^k)
  ("k-body" traceless operators).  dim V_k = C(N,k)^2 - C(N,k-1)^2.
* Adjoint Casimir C(X) = sum_{ij} [E_ij, [E_ji, X]], E_ij = psibar_i psi_j, has eigenvalue 2 k (N + 1 - k) on V_k
  (checked numerically in tests/test_size_decomposition.py).
* The normalized character of V_k at g in U(N) is chi_k(g) = (|chi_{Lambda^k}(g)|^2 - |chi_{Lambda^{k-1}}(g)|^2)/dim V_k.
  For the parity U_S = (-1)^{N_S} (|S| = s), chi_{Lambda^k}(U_S) = Krawtchouk K_k(s) = sum_j (-1)^j C(s,j) C(N-s,k-j).

Functions: casimir_eigs, dim_Vk, krawtchouk, chi_hat, size_weights (w_k = ||P^{(k)}||_F^2 / Tr P for Hermitian P).
"""
from __future__ import annotations

from math import comb

import numpy as np

from q_scan import subset_index


def dim_Vk(N, k):
    return comb(N, k) ** 2 - (comb(N, k - 1) ** 2 if k >= 1 else 0)


def casimir_eigs(N, P):
    return [2 * k * (N + 1 - k) for k in range(min(P, N - P) + 1)]


def krawtchouk(N, s, k):
    return sum((-1) ** j * comb(s, j) * comb(N - s, k - j) for j in range(0, min(s, k) + 1))


def chi_hat(N, s, k):
    """Normalized character of V_k at the parity U_S, |S| = s."""
    num = krawtchouk(N, s, k) ** 2 - (krawtchouk(N, s, k - 1) ** 2 if k >= 1 else 0)
    return num / dim_Vk(N, k)


class HopTable:
    """Index tables for E_ij = psibar_i psi_j on Lambda^P (bitmask basis of q_scan.subset_index)."""

    def __init__(self, N, P):
        self.N, self.P = N, P
        self.masks, index = subset_index(N, P)
        self.D = len(self.masks)
        self.ops = []
        for i in range(N):
            for j in range(N):
                src, dst, sgn = [], [], []
                for b, m in enumerate(self.masks):
                    if not (m >> j) & 1:
                        continue
                    if i == j:
                        src.append(b); dst.append(b); sgn.append(1.0)
                        continue
                    if (m >> i) & 1:
                        continue
                    s1 = bin(m & ((1 << j) - 1)).count("1")           # psi_j
                    m2 = m & ~(1 << j)
                    s2 = bin(m2 & ((1 << i) - 1)).count("1")          # psibar_i
                    src.append(b); dst.append(index[m2 | (1 << i)]); sgn.append((-1.0) ** (s1 + s2))
                self.ops.append((np.array(src), np.array(dst), np.array(sgn)))

    def M(self, X):
        """M(X) = sum_ij E_ij X E_ij^dagger."""
        out = np.zeros_like(X)
        for src, dst, sgn in self.ops:
            if len(src) == 0:
                continue
            out[np.ix_(dst, dst)] += (sgn[:, None] * sgn[None, :]) * X[np.ix_(src, src)]
        return out

    def casimir(self, X, cP=None):
        """C(X) = 2 (c_P X - M(X)), c_P = sum_ij E_ij E_ji on Lambda^P (a scalar, computed on X = 1 if not given)."""
        if cP is None:
            cP = self.cP()
        return 2 * (cP * X - self.M(X))

    def cP(self):
        # sum_ij E_ij E_ji = sum_ij E_ij E_ij^dagger restricted: its value on any state; use M(1)
        one = np.eye(self.D)
        return float(np.real(self.M(one)[0, 0]))


def size_weights(P_mat, table: HopTable):
    """w_k = ||Pi_k(P)||_F^2 / Tr(P) for Hermitian P on Lambda^P, via Lagrange filters in the adjoint Casimir.
    Returns (w, residual) with residual = |sum_k Pi_k(P) - P|_max (should be ~1e-10)."""
    N, Pp = table.N, table.P
    eigs = casimir_eigs(N, Pp)
    cP = table.cP()
    comps = []
    for k, ek in enumerate(eigs):
        V = P_mat.copy()
        for kk, ekk in enumerate(eigs):
            if kk == k:
                continue
            V = (table.casimir(V, cP) - ekk * V) / (ek - ekk)
        comps.append(V)
    trP = float(np.real(np.trace(P_mat)))
    w = np.array([float(np.real(np.vdot(V, V))) / trP for V in comps])
    resid = float(np.max(np.abs(sum(comps) - P_mat)))
    return w, resid


# ----------------------------------------------------------------------------- torus / character method
def overlap_class_averages(Pm, masks, N, P):
    """For a (Hermitian) operator Pm on Lambda^P in the bitmask basis, return, for r = 0..min(P, N-P):
       H[r] = mean |P_ST|^2      over ordered pairs with |S \\ T| = r   (off-diagonal content -> size weights w_k)
       G[r] = mean P_SS P_TT     over ordered pairs with |S \\ T| = r   (diagonal correlations -> transmissions phi_k)
    These are the S_N-orbit averages that determine the U(N)-twirls via characters on the maximal torus."""
    m = np.array(masks, dtype=np.int64)
    D = len(m)
    rmax = min(P, N - P)
    diag = np.real(np.diag(Pm))
    absq = np.abs(Pm) ** 2
    H = np.zeros(rmax + 1)
    G = np.zeros(rmax + 1)
    cnt = np.zeros(rmax + 1)
    # |S \ T| = P - |S cap T|; popcount of (m_S & m_T) row by row
    pc = np.vectorize(lambda v: bin(int(v)).count("1"))
    for s0 in range(0, D, 256):
        blk = m[s0:s0 + 256]
        inter = pc(blk[:, None] & m[None, :])
        r = P - inter
        for rr in range(rmax + 1):
            sel = r == rr
            H[rr] += absq[s0:s0 + 256][sel].sum()
            G[rr] += (diag[s0:s0 + 256, None] * diag[None, :])[sel].sum()
            cnt[rr] += sel.sum()
    return H / cnt, G / cnt, cnt


def character_coeffs(vals, N, P):
    """Solve sum_k c_k [C(N-2r,k-r) - C(N-2r,k-r-1)] = vals[r] C(N-2r,P-r) (r = 0..P) for the coefficients c_k of
    a class function sum_{S,T} f(|S\\T|) z^S zbar^T = sum_k c_k chi_{V_k}(z) (triangular in r <= k)."""
    rmax = len(vals) - 1
    A = np.zeros((rmax + 1, rmax + 1))
    for r in range(rmax + 1):
        for k in range(rmax + 1):
            t1 = comb(N - 2 * r, k - r) if k - r >= 0 else 0
            t2 = comb(N - 2 * r, k - r - 1) if k - r - 1 >= 0 else 0
            A[r, k] = t1 - t2
    rhs = np.array([vals[r] * comb(N - 2 * r, P - r) for r in range(rmax + 1)])
    return np.linalg.solve(A, rhs)


def size_and_transmission(Pm, masks, N, P):
    """Return (w_char, phi) for a projector Pm of rank d:
       w_char[k] = c^H_k dim V_k / d   (torus/S_N estimate of the size weights; = Casimir w_k in expectation)
       phi[k]    = c^G_k                (eigenvalue of X -> P X P on V_k, U(N)-twirled; phi_0 = a)"""
    d = float(np.real(np.trace(Pm)))
    H, G, _ = overlap_class_averages(Pm, masks, N, P)
    cH = character_coeffs(H, N, P)
    cG = character_coeffs(G, N, P)
    w = np.array([cH[k] * dim_Vk(N, k) / d for k in range(len(cH))])
    return w, cG


# ----------------------------------------------------------------------------- exact transmissions from size weights
def johnson_idempotent_values(N, P, masks=None):
    """e[k][r] = entry (E_k)_{ST} of the k-th primitive idempotent of the Johnson scheme J(N, P) at |S \\ T| = r.
    E_k projects functions on P-subsets onto the S_N-isotypic component (N-k, k) = the weight-zero part of V_k
    (the 'k-body' diagonal operators).  Computed numerically from the distance-1 adjacency."""
    if masks is None:
        masks, _ = subset_index(N, P)
    m = np.array(masks, dtype=np.int64)
    D = len(m)
    inter = np.array([[bin(int(x & y)).count("1") for y in m] for x in m])
    dist = P - inter
    A1 = (dist == 1).astype(float)
    ev, V = np.linalg.eigh(A1)
    rmax = min(P, N - P)
    e = np.zeros((rmax + 1, rmax + 1))
    for k in range(rmax + 1):
        theta = (P - k) * (N - P - k) - k
        sel = np.abs(ev - theta) < 1e-6
        Ek = V[:, sel] @ V[:, sel].T
        for r in range(rmax + 1):
            e[k, r] = Ek[dist == r].mean()
    return e


def transmission_matrix(N, P, masks=None):
    """Universal matrix A (rows j, cols k) with phi_j = sum_k A_jk c_k, c_k = ||P^{(k)}||^2 / dim V_k.
    Column k = character coefficients of the class function g -> ||rho(g)^{(k)}||^2 (torus: Johnson idempotent E_k)."""
    e = johnson_idempotent_values(N, P, masks)
    return np.column_stack([character_coeffs(e[k], N, P) for k in range(e.shape[0])])


def transmissions_from_weights(w, d, N, P, A=None):
    """Exact U(N)-twirled transmissions phi_j of X -> P X P (per realization) from the exact size weights w_k."""
    if A is None:
        A = transmission_matrix(N, P)
    c = np.array([w[k] * d / dim_Vk(N, k) for k in range(len(w))])
    return A @ c


# ----------------------------------------------------------------------------- decoder T_4: size-resolved remainder
def size_components(P_mat, table: HopTable):
    """Return the list of size components P^{(k)} (Lagrange filters in the adjoint Casimir) and the residual."""
    eigs = casimir_eigs(table.N, table.P)
    cP = table.cP()
    comps = []
    for k, ek in enumerate(eigs):
        V = P_mat.copy()
        for kk, ekk in enumerate(eigs):
            if kk != k:
                V = (table.casimir(V, cP) - ekk * V) / (ek - ekk)
        comps.append(V)
    return comps, float(np.max(np.abs(sum(comps) - P_mat)))


def decoder_T4_pieces(P_mat, comps, phi):
    """Exact per-realization pieces of T_4(U) = c_0 + 2a(T_2(U) - a) + R(U) (research note section 11):
       c_0   = sum_k phi_k w_k               (U(N)-twirled part = two-independent-model value)
       Gamma = [<P^{(k)}, delta_P P^{(k')}>]/d  for k, k' >= 1, delta_P = (X -> P X P) - twirl
    Returns (c0, Gamma, w)."""
    d = float(np.real(np.trace(P_mat)))
    w = np.array([float(np.real(np.vdot(C, C))) / d for C in comps])
    c0 = float(np.dot(phi, w))
    K = len(comps)
    G = np.zeros((K - 1, K - 1))
    for j in range(1, K):
        PCj = P_mat @ comps[j] @ P_mat                      # computed on the fly (memory)
        for i in range(1, K):
            val = float(np.real(np.vdot(comps[i], PCj))) / d
            if i == j:
                val -= phi[i] * w[i]
            G[i - 1, j - 1] = val
        del PCj
    return c0, G, w
