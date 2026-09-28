"""Tests for src/lmrs_predictions.py (derivation: docs/derivations.md D1).  Run: .venv/bin/python tests/test_lmrs_predictions.py"""
import os
import sys
from math import cos, gamma, pi

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from lmrs_predictions import ALPHA_S_Q3, b_q, r_charge, var_capped, var_zero_energy, zero_energy_expectation, zero_energy_factor  # noqa: E402


def test_reproduces_lmrs_table1():
    """LMRS Table 1 (N = 16, normalized by the total N_BPS over j = 0, +-1/3): psi 0.111, psibar psi 0.0874."""
    L = 1 + 2 * cos(pi / 3)
    pref = lambda D: (2 * ALPHA_S_Q3 * 16) ** (-2 * D)                                   # noqa: E731
    psi = b_q(3) * pref(1 / 6) * cos(pi / 3) / (pi * L) * gamma(0.5) * gamma(1 / 3 + 0.5)
    neutral = b_q(3) ** 2 * pref(1 / 3) * zero_energy_factor(1 / 3, 0) / L
    assert abs(psi - 0.111) < 5e-4 and abs(neutral - 0.0874) < 5e-5


def test_identity_normalization():
    """Delta -> 0 is the identity operator: tau_j(1) = 1 in every sector."""
    for j in (0.0, 1 / 6, -1 / 6, 1 / 3):
        assert abs(zero_energy_factor(1e-9, j) - 1) < 1e-6


def test_r_charge_assignment_and_parity_factor():
    assert r_charge(14, 7, 3) == 0 and abs(r_charge(15, 7, 3) + 1 / 6) < 1e-12
    ratio = zero_energy_factor(1 / 3, -1 / 6) / zero_energy_factor(1 / 3, 0)
    assert abs(ratio - 0.9204) < 1e-3


def test_wavefunction_integral_reproduces_gamma_formula():
    """<Z_j|e^{-Delta ell}|Z_j> from the LMRS eq. (22) wavefunctions equals the closed form of eq. (66)/Tr P_j."""
    import numpy as np
    for D in (1 / 6, 1 / 3, 2 / 3):
        for j in (0.0, -1 / 6, 1 / 3):
            assert abs(zero_energy_expectation(lambda l: np.exp(-D * l), j) / zero_energy_factor(D, j) - 1) < 1e-6


def test_cap_only_lowers_and_is_inactive_for_large_cap():
    assert var_capped(14, 0.25) < var_zero_energy(14)
    assert abs(var_capped(14, 1e6) / var_zero_energy(14) - 1) < 1e-6


def test_mellin_transform_of_K_squared():
    """GR 6.576.4 special case used in docs/derivations.md D1 step 4c."""
    import numpy as np
    from scipy.integrate import quad
    from scipy.special import gamma as G, kv
    for s, nu in [(2.0, 0.5), (8 / 3, 1 / 3), (8 / 3, 2 / 3), (3.0, 5 / 6), (2.2, 1 / 6)]:
        num = quad(lambda x: x ** (s - 1) * kv(nu, x) ** 2, 0, np.inf, limit=400)[0]
        cf = np.sqrt(np.pi) * G(s / 2) * G(s / 2 - nu) * G(s / 2 + nu) / (4 * G((s + 1) / 2))
        assert abs(num / cf - 1) < 1e-9


if __name__ == "__main__":
    for name, f in list(globals().items()):
        if name.startswith("test_"):
            f()
            print("ok", name)
