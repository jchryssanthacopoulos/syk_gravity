r"""
lmrs_predictions.py -- parameter-free N=2 super-Schwarzian predictions (Lin-Maldacena-Rozenberg-Shan 2022,
arXiv:2207.00408, "LMRS") for BPS-projected occupation operators in one-flavor N=2 SYK.  Derivation:
docs/derivations.md, section D1.

Quantity.  For the uplift decoder X = P_B Pi P_B, Pi = 1 - n_i, on the BPS space B at fixed fermion number P,
    Var(lambda) = tau(P n P n P) - tau(P n P)^2 = Tr[P_j O P_j O] / Tr[P_j],   O = n_i - mu,
i.e. exactly the zero-energy ("infinitely long time") two-point function of the neutral bilinear O in the single
R-charge sector j = (P - N/2)/q_hat that B occupies.  LMRS eqs. (66), (68), (85) give, with Delta = 2/(2 q_hat),

    Var_j = b_q^2 (2 alpha_S N)^{-2 Delta} * cos(pi j)/(2 pi) * Delta Gamma(Delta)^2 Gamma(Delta+1/2+j)
                                                                  Gamma(Delta+1/2-j) / Gamma(2 Delta)

    b_q = [tan(pi/(2 q_hat)) / (2 pi)]^{1/q_hat}     (LMRS eq. 78),   alpha_S = 0.00842 for q_hat = 3 (eq. 80).

Assumptions (see derivations.md): (i) large N: the same-index bilinear n_i - <n_i> has, at leading order in 1/N,
the same conformal two-point function as the i != k bilinear psi^k psibar^i treated in LMRS eq. (85) (one Wick
contraction, b_q^2/(J t)^{2 Delta}); (ii) Schwarzian regime (LMRS 'reliable' for the zero-energy sector at large N);
(iii) the Schwarzian coupling alpha_S is the large-N numerical value.  There are no fitted parameters.
Only q_hat = 3 has a known alpha_S; for other q_hat the function takes alpha_S as an argument.
"""
from __future__ import annotations

from math import cos, gamma, pi, tan

import numpy as np

ALPHA_S_Q3 = 0.00842          # LMRS eq. (80), q_hat = 3


def b_q(q: int) -> float:
    """Conformal fermion two-point normalization, LMRS eq. (78)."""
    return (tan(pi / (2 * q)) / (2 * pi)) ** (1.0 / q)


def r_charge(Np: int, P: int, q: int) -> float:
    """R-charge of the fermion-number-P sector of the N'-fermion theory: j = (P - N'/2)/q_hat.
    (Verified: BPS states occur only for |j| < 1/2, with dim ratio cos(pi j) -- see derivations.md D1.)"""
    return (P - Np / 2) / q


def zero_energy_factor(Delta: float, j: float) -> float:
    """cos(pi j)/(2 pi) * Delta Gamma(Delta)^2 Gamma(Delta+1/2+j) Gamma(Delta+1/2-j) / Gamma(2 Delta):
    LMRS eq. (66) divided by Tr P_j = e^{S0} cos(pi j) (eqs. 38, 68)."""
    return (cos(pi * j) / (2 * pi) * Delta * gamma(Delta) ** 2
            * gamma(Delta + 0.5 + j) * gamma(Delta + 0.5 - j) / gamma(2 * Delta))


def zero_energy_density(ell, j):
    """Normalized-up-to-constant density in the renormalized length ell of the zero-energy state |Z_j>
    (LMRS eq. 22): |g1|^2 + |g2|^2 ~ e^{-ell} [K_{1/2-j}^2 + K_{1/2+j}^2](2 e^{-ell/2})."""
    from scipy.special import kv
    z = 2 * np.exp(-ell / 2)
    return np.exp(-ell) * (kv(0.5 - j, z) ** 2 + kv(0.5 + j, z) ** 2)


def zero_energy_expectation(f, j, lo=-40.0, hi=400.0):
    """<Z_j| f(ell) |Z_j> / <Z_j|Z_j>.  With f = e^{-Delta ell} this reproduces zero_energy_factor(Delta, j).
    The density decays like e^{-(1/2 - |j|) ell} at large ell, hence the long upper limit."""
    from scipy.integrate import quad
    num = quad(lambda l: f(l) * zero_energy_density(l, j), lo, hi, limit=400)[0]
    den = quad(lambda l: zero_energy_density(l, j), lo, hi, limit=400)[0]
    return num / den


def var_capped(Np, cap, P=None, q=3, alpha_S=None):
    """TOY MODEL (docs/derivations.md D1, 'Why the mismatch'): the LMRS zero-energy average with the conformal
    matter correlator b_q^2 (J t)^{-2 Delta}, J t = 2 alpha_S N' e^{ell/2}, replaced by min(conformal, cap), where
    cap is the operator's exact UV equal-time connected value (b(1-b) for n_i).  Not a controlled approximation."""
    P = Np // 2 if P is None else P
    alpha_S = ALPHA_S_Q3 if alpha_S is None else alpha_S
    Delta, j = 2.0 / (2 * q), r_charge(Np, P, q)
    A = b_q(q) ** 2 * (2 * alpha_S * Np) ** (-2 * Delta)
    return zero_energy_expectation(lambda l: min(A * np.exp(-Delta * l), cap), j)


def var_zero_energy(Np: int, P: int | None = None, q: int = 3, alpha_S: float | None = None) -> float:
    """LMRS prediction for Var(lambda) of P_B (1 - n_i) P_B in the charge-P BPS sector of N'-fermion SYK."""
    P = Np // 2 if P is None else P
    if alpha_S is None:
        if q != 3:
            raise ValueError("alpha_S is only known for q_hat = 3; pass it explicitly")
        alpha_S = ALPHA_S_Q3
    Delta = 2.0 / (2 * q)
    j = r_charge(Np, P, q)
    if abs(j) >= 0.5:
        raise ValueError(f"no BPS states in the Schwarzian description for |j| = {abs(j)} >= 1/2")
    return b_q(q) ** 2 * (2 * alpha_S * Np) ** (-2 * Delta) * zero_energy_factor(Delta, j)


def r_pred(Np: int, P: int | None = None, q: int = 3, alpha_S: float | None = None) -> float:
    """Predicted r = Var(lambda) / [b(1-b)], b = 1 - P/N' (UV Bernoulli variance)."""
    P = Np // 2 if P is None else P
    b = 1 - P / Np
    return var_zero_energy(Np, P, q, alpha_S) / (b * (1 - b))


if __name__ == "__main__":
    # Reproduce LMRS Table 1 (N = 16, all j, normalized by the total N_BPS) as a check of the transcription.
    L = sum(cos(pi * j) for j in (0, 1 / 3, -1 / 3))
    pref = lambda D: (2 * ALPHA_S_Q3 * 16) ** (-2 * D)                                   # noqa: E731
    psi_j0 = b_q(3) * pref(1 / 6) * cos(0) * cos(-pi / 3) / (pi * L) * gamma(0.5) * gamma(1 / 3 + 0.5)
    neutral_j0 = b_q(3) ** 2 * pref(1 / 3) * zero_energy_factor(1 / 3, 0) * cos(0) / L
    print(f"LMRS Table 1 check: psi j=0 -> {psi_j0:.4f} (table 0.111); psibar psi j=0 -> {neutral_j0:.4f} (0.0874)")
    for Np in range(8, 21):
        print(Np, f"j={r_charge(Np, Np // 2, 3):+.3f}  Var={var_zero_energy(Np):.4f}  r={r_pred(Np):.4f}")
