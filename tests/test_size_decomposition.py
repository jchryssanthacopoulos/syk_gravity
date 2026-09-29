"""Tests for src/size_decomposition.py and src/parity_decoder.py (research/notes/decoder_moments_2026-09-29.md).
Run: .venv/bin/python tests/test_size_decomposition.py"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from chord_invariants import BPSProjector  # noqa: E402
from parity_decoder import parity_projector_diag, x_exact  # noqa: E402
from size_decomposition import HopTable, casimir_eigs, chi_hat, dim_Vk, size_weights  # noqa: E402


def _superop(T):
    D = T.D
    cP = T.cP()
    S = np.zeros((D * D, D * D))
    for c in range(D * D):
        X = np.zeros(D * D)
        X[c] = 1
        S[:, c] = T.casimir(X.reshape(D, D), cP).ravel()
    return S


def test_casimir_spectrum_and_characters():
    N, P = 6, 3
    T = HopTable(N, P)
    S = _superop(T)
    w, V = np.linalg.eig(S)
    w = w.real
    for k, e in enumerate(casimir_eigs(N, P)):
        assert int(np.sum(np.abs(w - e) < 1e-6)) == dim_Vk(N, k)
    Vinv = np.linalg.inv(V)
    for s in (1, 2):
        u = 2 * parity_projector_diag(T.masks, list(range(s))) - 1
        Ad = (u[:, None] * u[None, :]).ravel()
        for k, e in enumerate(casimir_eigs(N, P)):
            sel = np.abs(w - e) < 1e-6
            tr = np.real(np.trace((Vinv[sel] * Ad[None, :]) @ V[:, sel]))
            assert abs(tr / dim_Vk(N, k) - chi_hat(N, s, k)) < 1e-9


def test_w0_equals_a_and_weights_sum_to_one():
    for seed in (0, 1):
        pr = BPSProjector(8, 4, q=3, seed=seed, method="basis")
        w, resid = size_weights(pr.B @ pr.B.conj().T, HopTable(8, 4))
        assert resid < 1e-10 and abs(w.sum() - 1) < 1e-10 and abs(w[0] - pr.a) < 1e-12


def test_haar_sum_rule_reproduces_weingarten():
    """Haar: E w_k = (1-a) dim V_k/(D^2-1) (k>=1) => E T_2 = a + (1-a)(|Tr U|^2/D^2... ) ; at half filling
    E T_2 = (a D^2 - 1)/(D^2 - 1), the exact Weingarten value."""
    from math import comb
    N, P = 10, 5
    D = comb(N, P)
    a = 0.4
    pred = a + (1 - a) * sum(dim_Vk(N, k) * chi_hat(N, 1, k) for k in range(1, P + 1)) / (D * D - 1)
    assert abs(pred - (a * D * D - 1) / (D * D - 1)) < 1e-12


def test_x_zero_at_half():
    for N in (8, 10, 12, 14):
        assert abs(x_exact(N, 3, N // 2)) < 1e-15


if __name__ == "__main__":
    for name, f in list(globals().items()):
        if name.startswith("test_"):
            f()
            print("ok", name)
