"""Tests for src/dssyk_n2.py (BLY arXiv:2308.16283).  Run: .venv/bin/python tests/test_dssyk_n2.py"""
import os
import sys
from math import comb

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from dssyk_n2 import (bps_fraction_chord, bps_fraction_index, hh_length_distribution, r_chord,  # noqa: E402
                      two_point_chord, two_point_schwarzian_limit)
from lmrs_predictions import zero_energy_factor  # noqa: E402


def test_index_formula_matches_exact_bps_counts():
    """BLY (A.4) reproduces exact BPS dimensions of our model in every charge sector (index saturation)."""
    from q_scan import harmonic, rand_form
    for N in (7, 8, 9, 10, 11):
        C = rand_form(1, N, 3)
        for P in range(N + 1):
            if comb(N, P) > 500 or abs(P - N / 2) / 3 >= 0.5:
                continue
            d = harmonic(C, N, P, 3).shape[1]
            assert abs(bps_fraction_index(N, P) * comb(N, P) - d) < 1e-6, (N, P, d)


def test_two_point_normalization_and_limits():
    for j in (0.0, 1 / 6, 1 / 3):
        assert abs(two_point_chord(1.3, j, 1e-9) - 1) < 1e-6                 # Delta -> 0: identity operator
    for j in (0.0, 1 / 6):                                                      # lambda -> 0: BLY (4.9)
        lam = 1e-2
        ratio = two_point_chord(lam, j, 1 / 3) / two_point_schwarzian_limit(lam, j, 1 / 3)
        assert abs(ratio - 1) < 3e-2, ratio
        # (4.9) carries the LMRS j-factor
        assert abs(two_point_schwarzian_limit(0.1, j, 1 / 3) / (0.2 ** (2 / 3)) / zero_energy_factor(1 / 3, j) - 1) < 1e-12


def test_bps_fraction_chord_small_lambda():
    """(3.17) -> sqrt(pi/lam) 2 cos(pi j) e^{-pi^2/(4 lam)} (leading, Schwarzian) for small lambda."""
    from math import cos, exp, pi, sqrt
    lam, j = 0.5, 0.2
    lead = sqrt(pi / lam) * exp(lam * j * j) * 2 * cos(pi * j) * exp(-pi ** 2 / (4 * lam))
    assert abs(bps_fraction_chord(lam, j) / lead - 1) < 1e-6


def test_r_chord_values():
    assert abs(r_chord(14) - 0.6352) < 1e-3 and abs(r_chord(15) - 0.5648) < 1e-3


def test_free_compression_is_zero_length_term():
    """r = a + sum_{n>=1} P_n q^{2 Delta n}: probabilities sum to 1, P_0 = a, and the q^{2 Delta n} average is r."""
    import numpy as np
    from math import exp
    for N in (8, 12, 14, 16):
        lam = 18 / N
        P = hh_length_distribution(lam, 0.0)
        n = np.arange(len(P))
        assert abs(P.sum() - 1) < 1e-10
        assert abs(P[0] - bps_fraction_chord(lam, 0.0)) < 1e-14
        assert abs((P * exp(-lam) ** (2 * n / 3)).sum() - two_point_chord(lam, 0.0, 1 / 3)) < 1e-10


if __name__ == "__main__":
    for name, f in list(globals().items()):
        if name.startswith("test_"):
            f()
            print("ok", name)
