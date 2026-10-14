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

BPS -> continuum two-point function of a charged chiral single-chord operator (derived here, docs/derivations.md D7;
verified against the transfer-matrix engine src/chord_transfer.py, tests/test_chord_dressing.py).  The operator has
p' fermions, raises the charge by delta = p'/p = 2 Delta and carries qt = q^delta at the chord junction (BBNR sec. 7).
From the ground state of charge s_out (|s_out| < 1/2) it reaches charge s_in = s_out + delta.
chord_energy(lam, s, phi)                 Lambda_s(phi) = 2 q^{-1/2}[cosh lam(s - 1/2) - cos phi]  (B family, BBNR 5.3)
b_family_density(lam, s, phi)             (q^2, e^{+-2i phi}; q^2) q^{-s} / (2 pi Lambda_s)          (BBNR 5.11-5.12)
bps_continuum_kernel(lam, s_out, delta, phi)  A(phi) = <v(phi)| qt^N |Psi>/<0|Psi>  (squared matrix element / density)
bps_continuum_weight(lam, s_out, delta, phi)  f = rho_B A = Al-Salam-Chihara weight, a = q^{s_in-1/2}, b = q^{1/2+delta-s_out}
bps_bps_weight(lam, s_out, delta)         its point mass at zero energy when s_in < 1/2 (= BLY eq. 4.14)
time_integrated_two_point(lam, s_out, delta)  int f / Lambda (2phi1 closed form); x 1/kappa = int_0^inf G(tau) dtau
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


# -----------------------------------------------------------------------------------------------------------------------
# BPS -> continuum two-point function of a charged chiral single-chord operator (docs/derivations.md D7)
# -----------------------------------------------------------------------------------------------------------------------

def _log_poch_pm(a: float, Q: float, cosphi, tol: float = 1e-17):
    """log (a e^{i phi}, a e^{-i phi}; Q)_inf = sum_k log(1 - 2 a Q^k cos phi + a^2 Q^{2k}), vectorized in cos phi.
    Real for real a; requires a != 1 where cos phi = 1 (each factor is >= (1 - a Q^k)^2)."""
    c = np.asarray(cosphi, dtype=float)
    out = np.zeros_like(c)
    t = float(a)
    while abs(t) > tol:
        out += np.log1p(-2 * t * c + t * t)
        t *= Q
    return out


def chord_energy(lam: float, s: float, phi):
    """B-family (range of Q) chord energy at charge s, BBNR (5.3): 2 q^{-1/2} [cosh lam(s - 1/2) - cos phi].
    Note: BBNR's eigenvector coefficients (5.10), alpha_n ~ H_n(cos phi | q^2), belong to the eigenvalue with +cos phi;
    with (5.3) the eigenvector is alpha_n(pi - phi), i.e. H_n(-cos phi | q^2) (checked against (5.5) directly)."""
    q = exp(-lam)
    return 2 * q ** -0.5 * (np.cosh(lam * (s - 0.5)) - np.cos(phi))


def b_family_density(lam: float, s: float, phi):
    """Density of B-family chord eigenstates in phi in [0, pi] at charge s, normalized so that its integral is the
    B-family fraction of the charge sector: (q^2, e^{+-2i phi}; q^2) q^{-s} / (2 pi Lambda_s(phi)) (BBNR 5.11-5.12)."""
    q = exp(-lam)
    Q = q * q
    phi = np.asarray(phi, dtype=float)
    log_e2 = np.log(4 * np.sin(phi) ** 2) + _log_poch_pm(Q, Q, np.cos(2 * phi))     # (e^{+-2i phi}; Q)
    return np.exp(log_poch(Q, Q) + log_e2) * q ** (-s) / (2 * pi * chord_energy(lam, s, phi))


def bps_continuum_kernel(lam: float, s_out: float, delta: float, phi):
    """A(phi) = <v_{s_in}(phi)| qt^N |Psi_{s_out}> / <0|Psi_{s_out}>, qt = q^delta, s_in = s_out + delta, with
    <0|v(phi)> = 1:  A = (q^{2 delta}; q^2) / [(q^{s_in + 3/2} e^{+-i phi}, q^{1/2 + delta - s_out} e^{+-i phi}; q^2)].
    Derivation (D7): Q^dag Psi = 0 gives <b_m| qt^N |Psi> = (1 - q^{2 delta}) (qt/q)^m (q^2;q^2)_m x_m with the ground-state
    amplitudes x_m = q^{3m/2} H_m(-cosh lam(s_out + 1/2) | q^2)/(q^2;q^2)_m (BLY 3.6); then the q-Mehler formula."""
    q = exp(-lam)
    Q = q * q
    s_in = s_out + delta
    c = np.cos(np.asarray(phi, dtype=float))
    return np.exp(log_poch(q ** (2 * delta), Q) - _log_poch_pm(q ** (s_in + 1.5), Q, c)
                  - _log_poch_pm(q ** (0.5 + delta - s_out), Q, c))


def bps_continuum_weight(lam: float, s_out: float, delta: float, phi):
    """Spectral weight f(phi) of the normalized BPS -> charge-(s_out + delta) two-point function,
    G(tau) = int_0^pi f(phi) exp(-tau kappa Lambda_{s_in}(phi)) dphi  [+ the point mass bps_bps_weight at zero energy].
    f = b_family_density x bps_continuum_kernel = normalized Al-Salam-Chihara weight (base q^2, x = cos phi):
        f = (ab; q^2) (q^2, e^{+-2i phi}; q^2) / [2 pi (a e^{+-i phi}, b e^{+-i phi}; q^2)],
        a = q^{s_in - 1/2},  b = q^{1/2 + delta - s_out}  (ab = q^{2 delta}).
    Integrates to 1 when s_in > 1/2 (no BPS states at s_in)."""
    q = exp(-lam)
    Q = q * q
    s_in = s_out + delta
    phi = np.asarray(phi, dtype=float)
    a, b = q ** (s_in - 0.5), q ** (0.5 + delta - s_out)
    log_e2 = np.log(4 * np.sin(phi) ** 2) + _log_poch_pm(Q, Q, np.cos(2 * phi))
    return np.exp(log_poch(a * b, Q) + log_poch(Q, Q) + log_e2 - _log_poch_pm(a, Q, np.cos(phi))
                  - _log_poch_pm(b, Q, np.cos(phi))) / (2 * pi)


def bps_bps_weight(lam: float, s_out: float, delta: float) -> float:
    """Point mass of the Al-Salam-Chihara measure at zero energy (present iff s_in = s_out + delta < 1/2):
    (q^{1 - 2 s_in}; q^2) / (q^{1 - 2 s_out}; q^2).  Equals the normalized BPS -> BPS weight of BLY (4.14)."""
    s_in = s_out + delta
    if s_in >= 0.5:
        return 0.0
    q = exp(-lam)
    return exp(log_poch(q ** (1 - 2 * s_in), q * q) - log_poch(q ** (1 - 2 * s_out), q * q))


def time_integrated_two_point(lam: float, s_out: float, delta: float, tol: float = 1e-16) -> float:
    """int_0^pi f(phi) / Lambda_{s_in}(phi) dphi for s_in > 1/2 (multiply by 1/kappa for int_0^inf G(tau) dtau):
        q^{s_in} / (1 - q^{2 delta}) 2phi1(q^2, q^2; q^{2 + 2 delta}; q^2, q^{2 s_in - 1}).
    Stieltjes transform of the ASC measure at x_a = (a + 1/a)/2, where the ASC polynomials are (ab;q^2)_n / (2a)^n
    (Casoratian sum, D7)."""
    s_in = s_out + delta
    if s_in <= 0.5:
        raise ValueError("needs s_in > 1/2 (no zero-energy states in the target sector)")
    q = exp(-lam)
    Q = q * q
    z = q ** (2 * s_in - 1)
    tot, term, n = 0.0, 1.0, 0
    while term > tol * tot or n < 3:
        tot += term
        term *= (1 - Q ** (n + 1)) / (1 - q ** (2 + 2 * delta) * Q ** n) * z
        n += 1
    return q ** s_in / (1 - q ** (2 * delta)) * tot


def time_integrated_two_point_small_lambda(lam: float, s_out: float, delta: float) -> float:
    """Leading lambda -> 0 behaviour of time_integrated_two_point: pi / (2 sin(pi delta)) lam^{delta-2} (2 s_in - 1)^{delta-1}
    (Schwarzian region dominates; valid for 0 < delta < 1)."""
    s_in = s_out + delta
    return pi / (2 * sin(pi * delta)) * lam ** (delta - 2) * (2 * s_in - 1) ** (delta - 1)
