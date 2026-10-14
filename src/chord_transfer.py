r"""
chord_transfer.py -- transfer-matrix ("chord engine") implementation of the one-sided double-scaled N=2 chord Hilbert
space at fixed charge s (Berkooz, Brukner, Narovlansky, Raz arXiv:2003.04405 = BBNR, secs. 4.4-4.7, 5.1, 7.1).

Independent numerical check of the closed forms in src/dssyk_n2.py (docs/derivations.md D7).

States are alternating strings of chord ends, read top -> bottom:
    |0>,  |n,X> = XOXO..XO,  |n,O> = OXOX..OX  (bosonic, length 2n),
    |n+1/2,X> = XOX..X,  |n+1/2,O> = OXO..O     (fermionic, length 2n+1).
Q opens an O chord at the bottom with weight q^{+s} and closes X chords; Q^dag opens X with q^{-s} and closes O
(BBNR 4.16-4.17); non-alternating results are null (BBNR 4.7).  Inner product: BBNR (4.37).  T = Q Q^dag + Q^dag Q.
The engine reproduces BBNR (4.40)-(4.42), BLY (3.17) and BLY (4.7) (tests/test_chord_dressing.py).

Matter (BBNR sec. 7.1): a chiral single-chord operator with p' fermions shifts the charge inside the operator region by
delta = p'/p and inserts qt^N at the junction, qt = q^delta, N|n,.> = n|n,.>.
"""
from __future__ import annotations

import numpy as np


def qpoch(a: float, q: float, n: int) -> float:
    """Finite Pochhammer (a; q)_n."""
    out = 1.0
    for k in range(n):
        out *= 1 - a * q ** k
    return out


class ChordN2:
    """Chord Hilbert space truncated at string length lmax + 1 (so T restricted to length <= lmax is exact)."""

    def __init__(self, lam: float, s: float, lmax: int):
        self.lam, self.s, self.lmax = lam, s, lmax
        q = np.exp(-lam)
        self.q = q
        self.states = [()]
        for l in range(1, lmax + 2):
            for first in ("O", "X"):
                other = "X" if first == "O" else "O"
                self.states.append(tuple(first if k % 2 == 0 else other for k in range(l)))
        self.idx = {st: i for i, st in enumerate(self.states)}
        n = len(self.states)
        Q = np.zeros((n, n))
        Qd = np.zeros((n, n))
        for j, v in enumerate(self.states):
            for M, opench, closech, oth in ((Q, "O", "X", "O"), (Qd, "X", "O", "X")):
                w = v + (opench,)                                      # open a chord at the bottom
                if w in self.idx:
                    M[self.idx[w], j] += q ** (s if opench == "O" else -s)
                nclose = sum(1 for c in v if c == closech)             # close a chord of type closech at position i
                for i, c in enumerate(v):
                    if c != closech:
                        continue
                    below = len(v) - 1 - i
                    above_oth = sum(1 for cc in v[:i] if cc == oth)
                    w = v[:i] + v[i + 1:]
                    if w in self.idx:
                        M[self.idx[w], j] += (-1) ** below * q ** (-above_oth + (nclose - 1))
        self.Q, self.Qd = Q, Qd
        G = np.zeros((n, n))                                           # Gram matrix, BBNR (4.37)
        q2 = q * q
        for i, v in enumerate(self.states):
            l = len(v)
            if l == 0:
                G[i, i] = 1.0
            elif l % 2 == 0:
                m = l // 2
                c = qpoch(q2, q2, m - 1)
                G[i, i] = q ** (-m) * c
                G[i, self.idx[tuple("X" if ch == "O" else "O" for ch in v)]] = -c
            else:
                m = (l - 1) // 2
                G[i, i] = q ** ((-s if v[0] == "O" else s) - m) * qpoch(q2, q2, m)
        self.G = G
        self.T_full = Q @ Qd + Qd @ Q
        keep = [i for i, v in enumerate(self.states) if len(v) <= lmax]
        self.keep = np.array(keep)
        self.boson = np.array([i for i in keep if len(self.states[i]) % 2 == 0])
        self.N = np.array([len(v) / 2 for v in self.states])

    def adjoint_error(self) -> float:
        """max |G Q - Qd^T G| on lengths <= lmax (Qd is the G-adjoint of Q)."""
        k = self.keep
        return float(np.abs((self.G @ self.Q - self.Qd.T @ self.G)[np.ix_(k, k)]).max())

    def bosonic(self):
        """(T, Gram, N) restricted to the bosonic sector (lengths 0, 2, ..., <= lmax)."""
        b = self.boson
        return self.T_full[np.ix_(b, b)], self.G[np.ix_(b, b)], self.N[b]


def geig(T: np.ndarray, G: np.ndarray):
    """Eigen-decomposition of a G-self-adjoint T: eigenvalues and G-orthonormal eigenvectors (columns)."""
    R = np.linalg.cholesky(G).T
    Ri = np.linalg.inv(R)
    Tt = R @ T @ Ri
    w, U = np.linalg.eigh((Tt + Tt.T) / 2)
    return w, Ri @ U


def ground_state_amplitudes(q: float, s: float, M: int):
    """Zero mode Psi = |0> + sum_m (x_m |m,X> + o_m |m,O>) of T_s for |s| < 1/2, from Q Psi = Q^dag Psi = 0:
        q^s o_m + q^-1 o_{m+1} - q^m x_{m+1} = 0,   q^-s x_m + q^-1 x_{m+1} - q^m o_{m+1} = 0   (x_0 = o_0 = 1).
    Closed form (BLY 3.6-3.7): x_m = q^{3m/2} H_m(-cosh lam(s + 1/2) | q^2)/(q^2;q^2)_m (o_m: s -> -s)."""
    o = np.zeros(M + 1)
    x = np.zeros(M + 1)
    o[0] = x[0] = 1.0
    for m in range(M):
        pre = -q / (1 - q ** (2 * m + 2))
        o[m + 1] = pre * (q ** s * o[m] + q ** (m + 1 - s) * x[m])
        x[m + 1] = pre * (q ** (m + 1 + s) * o[m] + q ** (-s) * x[m])
    return x, o


def bps_matter_spectral_weights(lam: float, s_out: float, delta: float, lmax: int = 300):
    """Discrete (truncated) spectral decomposition of the normalized BPS -> charge (s_out + delta) two-point function of a
    chiral single-chord operator with qt = q^delta (BBNR 7.8 with the outer interval projected on the zero mode):
        G(u) = <0| e^{-u T_in} qt^N |Psi_out> / <0|Psi_out> = sum_k f_k exp(-u Lambda_k).
    Returns (Lambda_k, f_k).  Requires |s_out| < 1/2."""
    Co = ChordN2(lam, s_out, lmax)
    To, Go, No = Co.bosonic()
    wo, Vo = geig(To, Go)
    if abs(wo[0]) > 1e-8:
        raise ValueError(f"no zero mode at s_out = {s_out} (lowest eigenvalue {wo[0]})")
    Psi = Vo[:, 0]
    Ci = ChordN2(lam, s_out + delta, lmax)
    Ti, Gi, _ = Ci.bosonic()
    wi, Vi = geig(Ti, Gi)
    qt = Co.q ** delta
    f = (Gi @ Vi)[0, :] * (Vi.T @ Gi @ (qt ** No * Psi)) / (Go @ Psi)[0]
    return wi, f
