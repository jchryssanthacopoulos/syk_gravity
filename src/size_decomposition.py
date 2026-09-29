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
