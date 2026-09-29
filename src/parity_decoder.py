r"""
parity_decoder.py -- the decoder in "parity" form and its Z2 chord rule (research/notes/decoder_moments_2026-09-29.md).

Exact reformulation.  Pi = 1 - n_i = (1 + U)/2 with U = (-1)^{n_i} (unitary involution).  Hence
    P Pi P = (P + P U P)/2,   lambda = (1 + mu)/2,   mu in spec(P U P),
and the decoder moments are m_k = 2^{-k} sum_j C(k,j) T_j with T_j = tau((P U P)^j).

Generalization.  U_S = (-1)^{N_S}, N_S = sum_{i in S} n_i, |S| = s.  For a supercharge term psi_I (|I| = p),
U_S psi_I U_S = (-1)^{|I cap S|} psi_I.  In any trace with m insertions of U_S, moving the U_S together leaves
tr(U_S^m) times a sign (-1)^{|I cap S|} for every Q-chord whose endpoints are separated by an odd number of U_S.
Averaged over the chord's random index set:
    x_s = E[(-1)^{|I cap S|}] = sum_k (-1)^k C(s,k) C(N-s,p-k) / C(N,p)      (exact, finite N)
        ~ q^{s/p} = exp(-2 p s / N)                                         (double-scaled, s << N)
For odd p and |S| = N/2 (N even), x = 0 exactly (complement symmetry: U_{S^c} = (-1)^{N_F} U_S).
"""
from __future__ import annotations

from math import comb, exp

import numpy as np


def parity_projector_diag(masks, S):
    """0/1 diagonal of Pi_S = (1 + (-1)^{N_S})/2 (even parity on S) in the occupation basis given by bitmasks."""
    Smask = 0
    for i in S:
        Smask |= 1 << i
    return np.array([1.0 if bin(m & Smask).count("1") % 2 == 0 else 0.0 for m in masks])


def x_exact(N, p, s):
    return sum((-1) ** k * comb(s, k) * comb(N - s, p - k) for k in range(p + 1)) / comb(N, p)


def x_double_scaled(N, p, s):
    return exp(-2.0 * p * s / N)


def mu_moments_from_lambda(lam, K=8):
    mu = 2 * np.asarray(lam) - 1
    return [float(np.mean(mu ** m)) for m in range(1, K + 1)]
