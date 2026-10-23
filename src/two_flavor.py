r"""
two_flavor.py -- the symmetrized two-flavor N=2 SYK model (Chang-Chen-Sia-Yang arXiv:2412.06902 eq. 3.18; CV sec. 6.2)
and the machinery for its BPS uplift N -> N+1 (docs/derivations.md D8).

Model.  Complex fermions psi_i, chi_i (i = 0..N-1).  Convention of CCSY: psi_i, chi_i are *creation* operators in the
formulas of the paper; here we write c^dag explicitly.  The supercharge
    Q = sum_{i,j,k} C_ijk  psi_i^dag psi_j^dag chi_k ,      C_ijk symmetric in (j,k)      (CCSY 3.18, CV 82)
creates two psi's and annihilates one chi: charge (+2, -1) under (N_psi, N_chi); H = {Q, Q^dag} conserves both.
The symmetrization makes Q commute with the bilinear V = sum_i psi_i^dag chi_i^dag (CCSY 3.19), and the tower
V^n |Omega> (|Omega> = all empty) consists of monotone BPS states (CCSY 3.20-3.25).

Couplings: C = (G + G^{T_jk})/2 with G i.i.d. complex Gaussian, E|G|^2 = 2 (the symmetrization of CCSY 3.18 applied to
i.i.d. couplings).  Only ratios of couplings enter BPS spaces and dressing costs.

Modes.  Every fermion is a bit of a mask; mode positions are arbitrary (psi_pos[i], chi_pos[i]).  For an uplift
N -> N+1 we place the new pair LAST (positions 2N, 2N+1), so that the old model is the restriction to the lower 2N bits
and the occupation pattern T = (t_psi, t_chi) of the new pair is read off the top two bits.  Jordan-Wigner signs count
occupied modes below the acted-on mode.
Sectors: (a, b) = (N_psi, N_chi).
"""
from __future__ import annotations

import itertools
from math import comb

import numpy as np
import scipy.sparse as sp


def random_symmetric_couplings(N: int, rng) -> np.ndarray:
    """C_ijk = (G_ijk + G_ikj)/2, G i.i.d. complex Gaussian with E|G|^2 = 2."""
    G = rng.standard_normal((N, N, N)) + 1j * rng.standard_normal((N, N, N))
    return (G + G.transpose(0, 2, 1)) / 2


def random_generic_couplings(N: int, rng) -> np.ndarray:
    """Unsymmetrized (Heydeman-Turiaci-Zhao) couplings: i.i.d. complex Gaussian, E|C|^2 = 2."""
    return rng.standard_normal((N, N, N)) + 1j * rng.standard_normal((N, N, N))


def _popcount(x):
    x = np.asarray(x, dtype=np.int64)
    c = np.zeros_like(x)
    while np.any(x):
        c += x & 1
        x = x >> 1
    return c


class TwoFlavor:
    """Fock space of N_s psi- and N_s chi-modes at arbitrary bit positions, and sector blocks of operators."""

    def __init__(self, C: np.ndarray, psi_pos=None, chi_pos=None):
        self.C = np.asarray(C)
        Ns = self.C.shape[0]
        self.Ns = Ns
        self.psi_pos = list(range(Ns)) if psi_pos is None else list(psi_pos)
        self.chi_pos = list(range(Ns, 2 * Ns)) if chi_pos is None else list(chi_pos)
        self.M = max(self.psi_pos + self.chi_pos) + 1
        self._sectors = {}

    # ------------------------------------------------------------------ basis
    def sector(self, a: int, b: int) -> np.ndarray:
        """Sorted masks with N_psi = a, N_chi = b."""
        key = (a, b)
        if key not in self._sectors:
            if not (0 <= a <= self.Ns and 0 <= b <= self.Ns):
                self._sectors[key] = np.zeros(0, dtype=np.int64)
            else:
                pm = [sum(1 << self.psi_pos[i] for i in S) for S in itertools.combinations(range(self.Ns), a)]
                cm = [sum(1 << self.chi_pos[i] for i in S) for S in itertools.combinations(range(self.Ns), b)]
                self._sectors[key] = np.sort(np.array([p | c for p in pm for c in cm], dtype=np.int64))
        return self._sectors[key]

    def dim(self, a, b):
        return len(self.sector(a, b))

    # ------------------------------------------------------------------ operators
    @staticmethod
    def _apply(masks, ops):
        """Apply a product of fermion operators to basis masks.  ops = [(pos, dagger), ...] listed LEFT to RIGHT as in
        the operator product; the rightmost acts first.  Returns (ok, sign, newmasks)."""
        m = masks.copy()
        sign = np.ones(len(m))
        ok = np.ones(len(m), dtype=bool)
        for pos, dag in reversed(ops):
            bit = np.int64(1) << np.int64(pos)
            occ = (m & bit) != 0
            ok &= (~occ) if dag else occ
            below = _popcount(m & (bit - 1))
            sign *= np.where(below % 2 == 1, -1.0, 1.0)
            m = np.where(occ, m & ~bit, m | bit)
        return ok, sign, m

    def op_block(self, terms, src, dst):
        """Sparse matrix of sum_t coef_t * (operator string t) from sector src to sector dst.
        terms: list of (coef, ops) with ops as in _apply."""
        sm, dm = self.sector(*src), self.sector(*dst)
        rows, cols, vals = [], [], []
        if len(sm) and len(dm):
            idx = np.arange(len(sm))
            for coef, ops in terms:
                if coef == 0:
                    continue
                ok, sign, nm = self._apply(sm, ops)
                if not np.any(ok):
                    continue
                r = np.searchsorted(dm, nm[ok])
                good = (r < len(dm))
                good[good] &= dm[r[good]] == nm[ok][good]
                rows.append(r[good])
                cols.append(idx[ok][good])
                vals.append(coef * sign[ok][good])
        if rows:
            rows, cols, vals = np.concatenate(rows), np.concatenate(cols), np.concatenate(vals)
        else:
            rows = cols = np.zeros(0, dtype=int)
            vals = np.zeros(0)
        return sp.csr_matrix((vals, (rows, cols)), shape=(len(dm), len(sm)), dtype=complex)

    def Q_terms(self, sites=None):
        """Terms of Q = sum_{i<j,k} (C_ijk - C_jik) psi_i^dag psi_j^dag chi_k, optionally restricted to a site subset."""
        S = range(self.Ns) if sites is None else sites
        S = sorted(S)
        out = []
        for i, j in itertools.combinations(S, 2):
            for k in S:
                coef = self.C[i, j, k] - self.C[j, i, k]
                if coef != 0:
                    out.append((coef, [(self.psi_pos[i], True), (self.psi_pos[j], True), (self.chi_pos[k], False)]))
        return out

    def Q(self, a, b, sites=None):
        """Q: (a, b) -> (a+2, b-1)."""
        return self.op_block(self.Q_terms(sites), (a, b), (a + 2, b - 1))

    def V(self, a, b, sites=None):
        """V = sum_i psi_i^dag chi_i^dag: (a, b) -> (a+1, b+1)."""
        S = range(self.Ns) if sites is None else sites
        return self.op_block([(1.0, [(self.psi_pos[i], True), (self.chi_pos[i], True)]) for i in S],
                             (a, b), (a + 1, b + 1))

    # ------------------------------------------------------------------ BPS
    def bps(self, a, b, tol=1e-9, return_laplacian_spectrum=False):
        """Orthonormal basis (columns) of harmonic states in sector (a, b): ker Q_out cap ker Q_in^dag."""
        d = self.dim(a, b)
        if d == 0:
            return np.zeros((0, 0), complex)
        L = np.zeros((d, d), complex)
        if self.dim(a + 2, b - 1):
            Qo = self.Q(a, b).toarray()
            L += Qo.conj().T @ Qo
        if self.dim(a - 2, b + 1):
            Qi = self.Q(a - 2, b + 1).toarray()
            L += Qi @ Qi.conj().T
        w, U = np.linalg.eigh((L + L.conj().T) / 2)
        scale = max(1.0, w.max()) if len(w) else 1.0
        keep = w < tol * scale
        if return_laplacian_spectrum:
            return U[:, keep], w
        return U[:, keep]


def vacuum_tower(model: TwoFlavor, n: int, sites=None):
    """V^n |Omega> as a vector in sector (n, n) (unnormalized)."""
    v = np.ones(1, complex)                      # |Omega> in sector (0,0)
    for k in range(n):
        v = model.V(k, k, sites) @ v
    return v
