r"""
nullmodels.py -- Haar / free-compression (Wachter, MANOVA) null model utilities.

Wachter(t, b): law of the compression of Bern(b) by a free projection of trace t, normalized per compressed
(BPS) dimension.  Continuous part on [lam_-, lam_+], lam_pm = (sqrt(t(1-b)) +- sqrt(b(1-t)))^2, plus atoms
max(0, 1 - b/t) at 0 and max(0, 1 - (1-b)/t) at 1.

(Moved here from scripts/analyze_chord_invariants.py, which keeps its own copy for backward reproducibility.)
"""
from __future__ import annotations

from math import comb, pi, sqrt

import numpy as np


def wachter_edges(t, b):
    lm = (sqrt(t * (1 - b)) - sqrt(b * (1 - t))) ** 2
    lp = (sqrt(t * (1 - b)) + sqrt(b * (1 - t))) ** 2
    return lm, lp


def wachter_density(t, b, xs):
    lm, lp = wachter_edges(t, b)
    xs = np.asarray(xs, dtype=float)
    f = np.zeros_like(xs)
    m = (xs > lm) & (xs < lp)
    f[m] = np.sqrt((lp - xs[m]) * (xs[m] - lm)) / (2 * pi * t * xs[m] * (1 - xs[m]))
    return f


def wachter_cdf(t, b, xs):
    """CDF (per-BPS normalization, atoms included)."""
    lm, lp = wachter_edges(t, b)
    p0, p1 = max(0.0, 1 - b / t), max(0.0, 1 - (1 - b) / t)
    g = np.linspace(lm, lp, 20001)[1:-1]
    f = np.sqrt(np.clip((lp - g) * (g - lm), 0, None)) / (2 * pi * t * g * (1 - g))
    c = np.concatenate([[0.0], np.cumsum((f[1:] + f[:-1]) / 2 * np.diff(g))])
    c = c / c[-1] * (1 - p0 - p1)
    xs = np.asarray(xs, dtype=float)
    out = np.interp(xs, g, c, left=0.0, right=1 - p0 - p1) + p0 * (xs >= 0)
    return np.where(xs >= 1, 1.0, out)


def ks_wachter(lam, t, b, tol=1e-9):
    """KS distance between an empirical spectrum (atoms included) and Wachter(t, b)."""
    lam = np.sort(np.where(lam > 1 - tol, 1.0, np.where(lam < tol, 0.0, lam)))
    n = len(lam)
    F = wachter_cdf(t, b, lam)
    emp_hi = np.searchsorted(lam, lam, side="right") / n
    emp_lo = np.searchsorted(lam, lam, side="left") / n
    F_lo = wachter_cdf(t, b, lam - 1e-12)
    return float(max(np.max(np.abs(emp_hi - F)), np.max(np.abs(emp_lo - F_lo))))


def free_compression_mu_moments(a, b, K=8):
    """Raw moments E[mu^m], m = 1..K, of mu = 2 lambda - 1 for lambda ~ Wachter(a, b), i.e. of P U P with
    U = 2 Pi - 1 (Pi ~ Bern(b)) compressed by a free projection of trace a."""
    from chord_invariants import free_compression_moments
    lam_m = [1.0] + free_compression_moments(b, a, K)          # E[lam^k], k = 0..K
    out = []
    for m in range(1, K + 1):
        out.append(sum(comb(m, k) * 2 ** k * lam_m[k] * (-1) ** (m - k) for k in range(m + 1)))
    return out
