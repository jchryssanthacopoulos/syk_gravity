"""Tests for src/chord_chaos.py and src/plethysm.py (research/notes/chord_size_dictionary_2026-10-02.md).
Run: .venv/bin/python tests/test_chord_chaos.py"""
import os
import sys
from math import exp

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from chord_chaos import (EvenSizeProjector, Sector, TopSizeFilter, bps_basis_from_H,  # noqa: E402
                         bps_projector_hodge, moment_ratio_exact, pi1_ratio_exact, random_couplings)
from plethysm import max_rows  # noqa: E402
from size_decomposition import HopTable, size_weights  # noqa: E402


def test_hodge_projector_matches_eigh():
    for N in (8, 10):
        S = Sector(N, 3)
        c = random_couplings(np.random.default_rng(3), S.M)
        B, w, _ = bps_basis_from_H(S.hamiltonian(c))
        P, w2, d = bps_projector_hodge(S, c)
        assert d == B.shape[1] == 2 * 3 ** (N // 2 - 1)
        assert np.abs(B @ B.conj().T - P).max() < 1e-12
        assert np.abs(np.sort(w) - w2).max() < 1e-10 * w.max()


def test_exact_gaussian_constants_by_monte_carlo():
    S = Sector(8, 3)
    rng = np.random.default_rng(1)
    tr, nrm = [], []
    for _ in range(600):
        H = S.hamiltonian(random_couplings(rng, S.M))
        t = np.trace(H).real
        tr.append(t)
        nrm.append(np.linalg.norm(H - t / S.D * np.eye(S.D)) ** 2)
    assert abs(np.mean(tr) - S.E_trH()) < 4 * np.std(tr) / np.sqrt(len(tr))
    assert abs(np.mean(nrm) - S.E_normH_nonid()) < 4 * np.std(nrm) / np.sqrt(len(nrm))


def test_pi1_closed_form_rationals():
    vals = {8: (7, 27), 10: (24, 55), 12: (220, 351), 14: (52, 63)}
    for N, (u, v) in vals.items():
        r = pi1_ratio_exact(N, 3)
        assert (r.numerator, r.denominator) == (u, v)
    assert str(moment_ratio_exact(14, 3)) == "17/14"


def test_pi1_double_scaled_limit_is_BLY():
    """pi_1/pi_0 -> 2q/(1-q) (BLY eqs. 3.6-3.10 at n = 1) as p -> infinity at fixed lambda = 2p^2/N."""
    for p in (25, 35):
        N = round(2 * p * p / 1.0)
        N += N % 2
        lam = 2 * p * p / N
        q = exp(-lam)
        assert abs(float(pi1_ratio_exact(N, p)) / (2 * q / (1 - q)) - 1) < 1e-3


def test_chaos1_identity_and_H_sizes():
    """p = 3, half filling: H has sizes {0, 2} only, and <P, H^(2)> = -d Tr H / D per realization."""
    N = 8
    S = Sector(N, 3)
    T = HopTable(N, N // 2)
    c = random_couplings(np.random.default_rng(5), S.M)
    H = S.hamiltonian(c)
    fH, _ = size_weights(H, T)
    assert np.all(np.abs(fH[[1, 3, 4]]) < 1e-12 * fH.sum())
    P, _, d = bps_projector_hodge(S, c)
    trH = np.trace(H).real
    H2 = H - trH / S.D * np.eye(S.D)
    assert abs(np.vdot(H2, P).real + d * trH / S.D) < 1e-9 * abs(d * trH / S.D)


def test_even_size_projector_and_top_filter():
    N = 10
    S = Sector(N, 3)
    c = random_couplings(np.random.default_rng(9), S.M)
    P, _, d = bps_projector_hodge(S, c)
    proj = EvenSizeProjector(N, N // 2)
    comps, res = proj.components(P, check=True)
    assert res < 1e-12
    assert np.abs(sum(comps) - P).max() < 1e-12
    assert abs(np.vdot(comps[0], comps[0]).real / d - d / S.D) < 1e-12      # w_0 = a
    H = S.hamiltonian(c)
    H2 = H @ H
    X4 = TopSizeFilter(N, N // 2).component(H2, 4)
    c2, _ = proj.components(H2)
    assert np.abs(X4 - c2[proj.ks.index(4)]).max() < 1e-9 * np.abs(X4).max()


def test_plethysm_row_bound():
    for n, p in ((2, 3), (3, 3), (2, 5)):
        mr, _ = max_rows(n, p)
        assert mr == n * (p - 1) + 1


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
