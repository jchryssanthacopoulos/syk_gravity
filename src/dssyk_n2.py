r"""
dssyk_n2.py -- finite-lambda (double-scaled, super-chord) predictions for N=2 SYK from Boruch, Lin & Yan,
"Exploring supersymmetric wormholes in N=2 SYK with chords", arXiv:2308.16283 ("BLY").  Derivation/usage:
docs/derivations.md D2; literature/notes/boruch_lin_yan_2023.md.

Conventions.  p-body supercharge (our q_hat = p = 3), N complex fermions (our N'), double-scaling parameter
lambda = 2 p^2 / N, q = e^{-lambda} (BLY eq. 1.6).  R-charge j in units where Q has charge 1: j = (P - N/2)/p up to a
sign (BLY eq. 1.5 uses the opposite sign; all formulas here are even in j).  Pochhammer (a; q)_inf.

Functions
---------
bps_fraction_index(N, P, p)   exact finite-N BPS fraction d / C(N, P) from the refined index (BLY App. A, eq. A.4;
                              assumes index saturation -- verified against exact counts, tests/test_dssyk_n2.py)
bps_fraction_chord(lam, j)    double-scaled BPS fraction D(j) (BLY eq. 3.17)
two_point_chord(lam, j, Delta) normalized zero-temperature two-point function of a neutral operator of dimension
                              Delta in the charge-j ground-state sector, normalized to its infinite-temperature value
                              (BLY eqs. 4.7 / 3.11 = 4.8).  For O = n_i - mu this is r = Var(lambda)/[b(1-b)].
two_point_schwarzian_limit    its lambda -> 0 limit (BLY eq. 4.9): (2 lambda)^{2 Delta} x LMRS j-factor
"""
from __future__ import annotations

from math import comb, cos, exp, gamma, pi, sin

import numpy as np


def log_poch(a: float, q: float, tol: float = 1e-18) -> float:
    """log (a; q)_inf = sum_{k>=0} log(1 - a q^k), for 0 <= a < 1, 0 <= q < 1 (log space: no underflow as q -> 1)."""
    from math import log1p
    out, t = 0.0, a
    while t > tol:
        out += log1p(-t)
        t *= q
    return out


def poch(a: float, q: float) -> float:
    """(a; q)_inf."""
    return exp(log_poch(a, q))


def bps_fraction_index(N: int, P: int, p: int = 3) -> float:
    """Exact d/C(N,P) from the refined index (BLY A.4):
    D_hat(j) = (1/p) sum_{n=0}^{p-1} sin^N(pi n/p) e^{i pi j p} e^{-2 pi i j n},  j = (P - N/2)/p,  d = D_hat 2^N.
    The inverse DFT only resolves j modulo 1, so this is valid for the representative sector |j| < 1/2 (where all
    BPS states live once the index is saturated); other sectors raise ValueError."""
    j = (P - N / 2) / p
    if abs(j) >= 0.5:
        raise ValueError(f"(A.4) applies only to |j| < 1/2 (got j = {j})")
    Dh = sum(sin(pi * n / p) ** N * np.exp(1j * pi * j * p) * np.exp(-2j * pi * j * n) for n in range(p)).real / p
    return float(Dh * 2 ** N / comb(N, P))


def bps_fraction_chord(lam: float, j: float) -> float:
    """BLY (3.17): D(j) = (q^{1+2j}; q^2)(q^{1-2j}; q^2)(q^2; q^2), q = e^{-lambda}."""
    q = exp(-lam)
    q2 = q * q
    return exp(log_poch(q ** (1 + 2 * j), q2) + log_poch(q ** (1 - 2 * j), q2) + log_poch(q2, q2))


def two_point_chord(lam: float, j: float, Delta: float) -> float:
    """BLY (4.7)/(3.11): <Psi_j| q^{2 Delta n} |Psi_j> / <Psi_j|Psi_j>
       = (q^{2+4D}; q^2) (q^{1+2j}; q^2)(q^{1-2j}; q^2)(q^2; q^2) / [ (q^{2+2D}; q^2)^2 (q^{1+2D+2j}; q^2)(q^{1+2D-2j}; q^2) ].
    Normalized so that the infinite-temperature (n = 0) value is 1."""
    q = exp(-lam)
    q2 = q * q
    L = lambda a: log_poch(a, q2)                                                          # noqa: E731
    return exp(L(q ** (2 + 4 * Delta)) + L(q ** (1 + 2 * j)) + L(q ** (1 - 2 * j)) + L(q2)
               - 2 * L(q ** (2 + 2 * Delta)) - L(q ** (1 + 2 * Delta + 2 * j)) - L(q ** (1 + 2 * Delta - 2 * j)))


def two_point_schwarzian_limit(lam: float, j: float, Delta: float) -> float:
    """BLY (4.9): (2 lambda)^{2 Delta} cos(pi j)/(2 pi) Delta Gamma(Delta)^2 Gamma(Delta + 1/2 +- j) / Gamma(2 Delta)."""
    return (2 * lam) ** (2 * Delta) * cos(pi * j) / (2 * pi) * Delta * gamma(Delta) ** 2 \
        * gamma(Delta + 0.5 + j) * gamma(Delta + 0.5 - j) / gamma(2 * Delta)


def r_chord(Np: int, P: int | None = None, p: int = 3) -> float:
    """Prediction for r = Var(lambda)/[b(1-b)] of the single-mode decoder P_B (1 - n_i) P_B at size N', charge P:
    lambda = 2 p^2 / N', Delta = 1/p (bilinear psibar_i psi_i: p_O = p_X = 1, BLY eq. 1.8), j = (P - N'/2)/p."""
    P = Np // 2 if P is None else P
    return two_point_chord(2 * p * p / Np, (P - Np / 2) / p, 1.0 / p)


def hh_length_distribution(lam: float, j: float, nmax: int = 2000, tol: float = 1e-18):
    """Chord-number (wormhole-length) distribution P_n of the supersymmetric Hartle-Hawking state |Psi, j>
    (BLY eqs. 3.2-3.3 recursion with alpha_0 = beta_0 = 1, weights 3.10, norm 3.11).  P_0 = D(j) is the weight of the
    empty wormhole |Omega> (eq. 3.15); P_n for n >= 1 is the n-th term of (3.10) divided by the norm.
    Terms are evaluated with rescaled amplitudes a q^{-n/2}, b q^{-n/2} (no overflow) and the series is truncated once
    terms fall below tol.  Returns numpy array P[0..n_last]."""
    q = exp(-lam)
    q2 = q * q
    a, b = 1.0, 1.0
    poch = 1.0                                          # (q^2; q^2)_{n-1}, updated incrementally
    terms = [0.0]
    for n in range(1, nmax):
        M = np.array([[-q ** (n - 1), 1 / q], [1 / q, -q ** (n - 1)]])
        a, b = np.linalg.solve(M, np.array([-q ** j * b, -q ** (-j) * a]))
        if n >= 2:
            poch *= 1 - q2 ** (n - 1)
        at, bt = a * q ** (-n / 2), b * q ** (-n / 2)
        t = poch * (at * at + bt * bt - 2 * at * bt * q ** n)
        terms.append(t)
        if n > 5 and abs(t) < tol:
            break
    norm = exp(-(log_poch(q ** (1 + 2 * j), q2) + log_poch(q ** (1 - 2 * j), q2) + log_poch(q2, q2)))
    P = np.array(terms) / norm
    P[0] = bps_fraction_chord(lam, j)
    return P
