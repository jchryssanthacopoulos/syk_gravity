r"""
chord_chaos.py -- exact finite-N chord number of the BPS projector (Wiener-chaos degree in the Gaussian couplings)
and its joint distribution with U(N) operator size.  Derivation: research/notes/chord_size_dictionary_2026-10-02.md
(summary: docs/derivations.md, D4).

Setting.  One-flavor N=2 SYK on N complex modes, supercharge Q = omega ^ (.) with omega a complex p-form whose
coefficients C_I (I a p-subset) are i.i.d. complex Gaussians, E|C_I|^2 = 2 (real and imaginary parts N(0,1); the
convention of q_scan.rand_form).  H = {Q, Q^dag} on the charge sector Lambda^F; P = BPS projector (kernel of H).

Chord number.  P(omega) is a square-integrable End(Lambda^F)-valued function of the Gaussian vector omega, so it has a
Wiener-chaos (Wick) expansion P = sum_n P_[n], P_[n] of bidegree (n, n) in (omega, omega-bar).  In the chord expansion
of E tau(P(omega) P(omega')) for correlated couplings E[C_I conj(C'_J)] = 2 x delta_IJ, P_[n] collects the diagrams
with n Q-chords and n Q^dag-chords crossing the cut, each weighted by x:
        E tau(P(omega) P(omega_x)) = sum_n pi_n x^(2n),     pi_n = E||P_[n]||_F^2 / d          (Mehler formula)
with omega_x = x omega + sqrt(1 - x^2) omega'.  pi_n is the exact finite-N analogue of BLY's chord-number distribution
P_n (arXiv:2308.16283, eqs. 3.10, 4.5-4.6).  The chaos grading commutes with U(N), so P has a joint (chord n, size k)
decomposition J_nk = E||P_[n]^(k)||^2 / d with marginals pi_n and w_k.

Exact results used here (p = 3, half filling):  P_[0] = a 1 (size 0);  P_[1] = gamma H^(2), with
gamma = -d (E Tr H / D) / E||H^(2)||^2 (Schur: the chaos-1, size-2 invariant space is one-dimensional), hence
pi_1 = d (E Tr H / D)^2 / E||H^(2)||^2 in closed form (pi1_closed).

Numerics.  Angle trick: for two independent coupling vectors omega, omega', the K couplings
omega_j = cos(theta_j) omega + sin(theta_j) omega', theta_j = j pi / K, are pairwise OU-correlated with
x = cos(theta_i - theta_j).  K diagonalizations give K^2 overlaps.  The size-resolved overlaps are estimated with
the U(N)/torus class-average method of size_decomposition.character_coeffs (unbiased in expectation).
"""
from __future__ import annotations

import itertools
from math import comb

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp

from q_scan import merge_sign, subset_index


# ----------------------------------------------------------------------------- couplings and Hamiltonian
def p_subsets(N, p):
    """Ordered list of p-subsets (tuples) and their bitmasks; coupling vectors are indexed in this order."""
    subs = list(itertools.combinations(range(N), p))
    masks = [sum(1 << i for i in T) for T in subs]
    return subs, masks


def random_couplings(rng, M):
    """Complex Gaussian coupling vector, E|C|^2 = 2 (same law as q_scan.rand_form)."""
    return rng.standard_normal(M) + 1j * rng.standard_normal(M)


class WedgeStructure:
    """Sparse structure of omega ^ (.) : Lambda^k -> Lambda^{k+p}; values are filled per coupling vector."""

    def __init__(self, N, k, p):
        self.N, self.k, self.p = N, k, p
        self.shape = (comb(N, k + p) if 0 <= k + p <= N else 0, comb(N, k) if 0 <= k <= N else 0)
        rows, cols, tidx, sgn = [], [], [], []
        if self.shape[0] and self.shape[1]:
            cmasks, _ = subset_index(N, k)
            _, rmap = subset_index(N, k + p)
            _, tmasks = p_subsets(N, p)
            for cj, S in enumerate(cmasks):
                for t, T in enumerate(tmasks):
                    if S & T:
                        continue
                    rows.append(rmap[S | T]); cols.append(cj); tidx.append(t); sgn.append(merge_sign(T, S))
        self.rows = np.array(rows, dtype=np.int64)
        self.cols = np.array(cols, dtype=np.int64)
        self.tidx = np.array(tidx, dtype=np.int64)
        self.sgn = np.array(sgn, dtype=float)

    def matrix(self, c):
        return sp.csr_matrix((self.sgn * c[self.tidx], (self.rows, self.cols)), shape=self.shape, dtype=complex)


def moment_ratio_exact(N, p, F=None):
    """m_2/m_1^2 - 1 (Fraction) for H = {Q, Q^dag} on Lambda^F, m_k = E tr H^k (normalized trace), p odd.
    Uses Tr_{Lambda^F} H = c |omega|^2 and sum_{IJ} ||{psi_I, psi_J^dag}||^2 = Sigma (see Sector.sigma_sum)."""
    from fractions import Fraction
    F = N // 2 if F is None else F
    D, M = comb(N, F), comb(N, p)
    h0 = comb(F, p) + comb(N - F, p)                      # scalar value of sum_I {psi_I, psi_I^dag} on Lambda^F
    Sig = sum(M * comb(p, m) * comb(N - p, p - m) * (comb(N - 2 * p + m, F - p + m) + comb(N - 2 * p + m, F - p))
              for m in range(1, p + 1))
    return Fraction(Sig, D * h0 * h0)


def pi1_ratio_exact(N, p, F=None):
    """pi_1 / pi_0 = m_1^2 / (m_2 - (1 + 1/M) m_1^2) as a Fraction.  Exact chord-1 weight for p = 3 at half filling
    (P_[1] = gamma H^(2) by Schur); for other p it is the weight of the projection of P_[1] onto :H:."""
    from fractions import Fraction
    return 1 / (moment_ratio_exact(N, p, F) - Fraction(1, comb(N, p)))

class Sector:
    """Charge-F sector of the p-form model on N modes: Hamiltonian builder, exact constants, size-class index."""

    def __init__(self, N, p, F=None):
        self.N, self.p = N, p
        self.F = N // 2 if F is None else F
        self.masks, _ = subset_index(N, self.F)
        self.D = len(self.masks)
        self.M = comb(N, p)
        self.up = WedgeStructure(N, self.F, p)            # Q on Lambda^F
        self.down = WedgeStructure(N, self.F - p, p)      # Q on Lambda^{F-p}; its adjoint maps into Lambda^F
        self.rmax = min(self.F, N - self.F)
        self._R = None

    def hamiltonian(self, c):
        """Dense H = Q^dag Q + Q Q^dag on Lambda^F for coupling vector c (length M)."""
        Qu = self.up.matrix(c)
        H = (Qu.conj().T @ Qu)
        if self.down.shape[0]:
            Qd = self.down.matrix(c)
            H = H + Qd @ Qd.conj().T
        return H.toarray()

    # exact Gaussian constants (E|C|^2 = 2) -------------------------------------------------------------------
    def c_diag(self):
        """Tr_{Lambda^F} {psi_I, psi_I^dag} (independent of I)."""
        N, p, F = self.N, self.p, self.F
        return comb(N - p, F) + comb(N - p, F - p)

    def E_trH(self):
        return 2.0 * self.M * self.c_diag()

    def sigma_sum(self):
        """Sum over ordered pairs (I, J) of ||{psi_I, psi_J^dag}||_F^2 on Lambda^F (odd p: zero unless I, J overlap)."""
        N, p, F, M = self.N, self.p, self.F, self.M
        tot = 0
        for m in range(1, p + 1):
            npairs = M * comb(p, m) * comb(N - p, p - m)
            tot += npairs * (comb(N - 2 * p + m, F - p + m) + comb(N - 2 * p + m, F - p))
        return tot

    def E_normH_nonid(self):
        """E ||H - (Tr H / D) 1||_F^2 (= E||H^(2)||^2 for p = 3 at half filling)."""
        c = self.c_diag()
        return 4.0 * (self.sigma_sum() - c * c * self.M / self.D)

    def pi1_closed(self, d):
        """Exact (pi_1, gamma) for p = 3 at half filling (P_[1] = gamma H^(2)); for other p, the share of the
        chaos-1 weight along the H direction (a lower bound on pi_1)."""
        EtrH = self.E_trH()
        EnH = self.E_normH_nonid()
        gamma = -d * (EtrH / self.D) / EnH
        return gamma * gamma * EnH / d, gamma

    # size classes ---------------------------------------------------------------------------------------------
    @property
    def R(self):
        """R[S, T] = |S \\ T| for basis subsets S, T of Lambda^F (int8)."""
        if self._R is None:
            m = np.array(self.masks, dtype=np.int64)
            pc = np.zeros((self.D, self.D), dtype=np.int8)
            for s0 in range(0, self.D, 512):
                inter = m[s0:s0 + 512, None] & m[None, :]
                pc[s0:s0 + 512] = np.vectorize(int.bit_count, otypes=[np.int64])(inter)
            self._R = (self.F - pc).astype(np.int8)
            self._Rflat = self._R.ravel()
            self._Rcount = np.bincount(self._Rflat, minlength=self.rmax + 1).astype(float)
        return self._R

    def class_sums(self, X, Y):
        """S[r] = sum over (S, T) with |S\\T| = r of Re(conj(X_ST) Y_ST)."""
        _ = self.R
        v = (np.conj(X) * Y).real.ravel()
        return np.bincount(self._Rflat, weights=v, minlength=self.rmax + 1)

    def class_counts(self):
        _ = self.R
        return self._Rcount


def bps_projector_hodge(S, c, d_expected=None, rel_tol=1e-9):
    """BPS projector on Lambda^F from the Hodge decomposition Lambda^F = ker H (+) im Q^dag (+) im Q (orthogonal since
    Q^2 = 0):  P = 1 - R1 R1^dag - R2 R2^dag, with R1 an orthonormal basis of the row space of Q_F (= im Q_F^dag) and
    R2 of the column space of Q_{F-p} (= im Q_{F-p}).  Two thin SVDs instead of a full eigendecomposition of H.
    Returns (P dense, spectrum of H on Lambda^F (ascending), d)."""
    D = S.D
    Qu = S.up.matrix(c).toarray()                         # (C(N,F+p) x D)
    _, s1, Vh1 = np.linalg.svd(Qu, full_matrices=False)
    del Qu
    tol = rel_tol * max(s1.max() if len(s1) else 0.0, 1.0)
    r1 = int((s1 > tol).sum())
    R = [Vh1[:r1].conj().T]
    spec = [s1[:r1] ** 2]
    if S.down.shape[0]:
        Qd = S.down.matrix(c).toarray()                   # (D x C(N,F-p))
        U2, s2, _ = np.linalg.svd(Qd, full_matrices=False)
        del Qd
        tol = rel_tol * max(s2.max(), s1.max() if len(s1) else 0.0, 1.0)
        r2 = int((s2 > tol).sum())
        R.append(U2[:, :r2])
        spec.append(s2[:r2] ** 2)
    d = D - sum(x.shape[1] for x in R)
    if d_expected is not None and d != d_expected:
        raise RuntimeError(f"BPS dimension {d} != expected {d_expected}")
    P = np.eye(D, dtype=complex)
    for Rx in R:
        P -= Rx @ Rx.conj().T
    w = np.sort(np.concatenate([np.zeros(d)] + spec))
    return P, w, d


def bps_basis_from_H(H, d_expected=None, rel_tol=1e-9):
    """Orthonormal basis of ker H (dense Hermitian); verifies a clean gap and, if given, the expected dimension."""
    w, Z = sla.eigh(H, overwrite_a=False, check_finite=False)
    scale = max(abs(w[-1]), 1.0)
    d = int((w < rel_tol * scale).sum())
    if d_expected is not None and d != d_expected:
        raise RuntimeError(f"BPS dimension {d} != expected {d_expected} (eigs near 0: {w[max(d-2,0):d+2]})")
    gap_ratio = (w[d] / scale) / max(abs(w[d - 1]) / scale, 1e-300) if 0 < d < len(w) else np.inf
    return Z[:, :d], w, gap_ratio


class EvenSizeProjector:
    """Exact U(N) size components of an operator on Lambda^F that has only EVEN sizes (true for the BPS projector and
    for H at half filling, p odd): with C~ = C / C_max, X_m = C~^m X (m = 0..r-1, r = number of even sizes) and the
    Vandermonde V_mk = (lambda_k / C_max)^m, X^(k) = sum_m (V^-1)_km X_m.  check=True applies C~ once more and returns
    the residual of the even-size assumption."""

    def __init__(self, N, F):
        from size_decomposition import HopTable, casimir_eigs
        self.table = HopTable(N, F)
        self.cP = self.table.cP()
        lam = np.array(casimir_eigs(N, F), dtype=float)
        self.ks = list(range(0, len(lam), 2))
        self.cmax = lam.max()
        nodes = lam[self.ks] / self.cmax
        r = len(self.ks)
        self.nodes = nodes
        self.Vinv = np.linalg.inv(np.vander(nodes, r, increasing=True).T)   # (V^-1)_{k m}

    def components(self, X, check=False):
        Xs = [X]
        r = len(self.ks)
        for _ in range(r - 1 + (1 if check else 0)):
            Xs.append(self.table.casimir(Xs[-1], self.cP) / self.cmax)
        comps = [sum(self.Vinv[i, m] * Xs[m] for m in range(r)) for i in range(r)]
        resid = None
        if check:
            pred = sum((self.nodes[i] ** r) * comps[i] for i in range(r))
            resid = float(np.max(np.abs(Xs[r] - pred)))
        return comps, resid


class TopSizeFilter:
    """Size-k component of an operator on Lambda^F whose sizes are EVEN and <= k (e.g. H^n at half filling, p = 3,
    k = 2n): product-form Lagrange filter prod_{k' < k, k' even} (C - lambda_k') / (lambda_k - lambda_k').  Uses k/2
    Casimir applications and O(1) dense matrices."""

    def __init__(self, N, F):
        from size_decomposition import HopTable, casimir_eigs
        self.table = HopTable(N, F)
        self.cP = self.table.cP()
        self.lam = casimir_eigs(N, F)

    def component(self, X, k):
        for kk in range(0, k, 2):
            X = (self.table.casimir(X, self.cP) - self.lam[kk] * X) / (self.lam[k] - self.lam[kk])
        return X

def size_overlaps_from_class_means(means, N, F, dimV=True):
    """Convert class means m[r] (of Re conj(X_ST) Y_ST) into E<X^(k), Y^(k)> = c_k dim V_k."""
    from size_decomposition import character_coeffs, dim_Vk
    c = character_coeffs(np.asarray(means), N, F)
    return np.array([c[k] * (dim_Vk(N, k) if dimV else 1) for k in range(len(c))])
