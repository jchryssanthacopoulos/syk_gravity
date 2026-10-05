"""Tests for scripts/run_ou_fourarc.py and scripts/analyze_ou_fourarc.py (docs/derivations.md D5)."""
import os
import sys

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))

from chord_chaos import Sector, bps_projector_hodge, p_subsets, random_couplings  # noqa: E402
from run_ou_fourarc import bps_basis, bps_dim  # noqa: E402
from analyze_ou_fourarc import krawtchouk_x, liberation_curve  # noqa: E402


def test_bps_dim_matches_known():
    # p = 3 half filling: 2 * 3^(N/2 - 1); odd N: 3^((N-1)/2)
    for N in (6, 8, 10, 12, 14):
        assert bps_dim(N, 3, N // 2) == 2 * 3 ** (N // 2 - 1)
    for N in (7, 9, 11, 13):
        assert bps_dim(N, 3, N // 2) == 3 ** ((N - 1) // 2)


def test_parity_conjugation_equals_sign_flipped_couplings():
    """D5 Sec. 1: U_A P(C) U_A = P(sigma C) with sigma_J = (-1)^{|A cap J|} -- the decoder is a coupling family."""
    N, p = 8, 3
    S = Sector(N, p)
    rng = np.random.default_rng(3)
    c = random_couplings(rng, S.M)
    _, tm = p_subsets(N, p)
    A = (1 << 2) | (1 << 5)
    sigma = np.array([(-1) ** bin(t & A).count("1") for t in tm], dtype=float)
    P, _, _ = bps_projector_hodge(S, c)
    Pf, _, _ = bps_projector_hodge(S, sigma * c)
    u = np.array([(-1) ** bin(m & A).count("1") for m in S.masks], dtype=float)
    assert np.abs(u[:, None] * P * u[None, :] - Pf).max() < 1e-10


def test_basis_and_moment_identity():
    N, p = 8, 3
    S = Sector(N, p)
    rng = np.random.default_rng(5)
    d = bps_dim(N, p, N // 2)
    V1, e1, _ = bps_basis(S, random_couplings(rng, S.M), d, rng)
    V2, e2, _ = bps_basis(S, random_couplings(rng, S.M), d, rng)
    assert max(e1, e2) < 1e-12 and np.abs(V1.conj().T @ V1 - np.eye(d)).max() < 1e-12
    c2 = np.linalg.svd(V1.conj().T @ V2, compute_uv=False) ** 2
    P1, P2 = V1 @ V1.conj().T, V2 @ V2.conj().T
    T2 = np.trace(P1 @ P2).real / d
    T4 = np.trace(P1 @ P2 @ P1 @ P2).real / d
    assert abs(T2 - c2.mean()) < 1e-12 and abs(T4 - (c2 ** 2).mean()) < 1e-12
    # T4 - 2 T2 + 1 = mean sin^4 >= 0
    assert abs((T4 - 2 * T2 + 1) - ((1 - c2) ** 2).mean()) < 1e-12


def test_liberation_closed_form():
    a = 0.4248
    t, m = liberation_curve(a, 2)
    u, b = np.exp(-t), 1 - a
    E = 1 - a * a * (2 - a) - 4 * a * b * b
    assert np.abs(m[0] - (a + b * u)).max() < 1e-9
    assert np.abs(m[1] - (a * a * (2 - a) + 4 * a * b * b * u + (E - 2 * a * b * b * t) * u * u)).max() < 1e-9


def test_krawtchouk_single_site():
    for N in (8, 12, 14):
        assert abs(krawtchouk_x(N, 3, 1) - (1 - 6 / N)) < 1e-14


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
