"""Tests for the chord dictionary of the two-flavor model and its uplift (src/two_flavor_chords.py; docs/derivations.md D9).
Plain-script runner: .venv/bin/python tests/test_two_flavor_chords.py"""
import os
import sys

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from two_flavor import TwoFlavor, random_generic_couplings, random_symmetric_couplings  # noqa: E402
from two_flavor_chords import DELTA, bps_fraction, dressing_targets, lam, resonant_terms, s_eff, two_point  # noqa: E402


def laplacian_spectrum(M, a, b):
    d = M.dim(a, b)
    L = np.zeros((d, d), complex)
    wB = np.zeros(0)
    if M.dim(a + 2, b - 1):
        Qo = M.Q(a, b).toarray()
        L += Qo.conj().T @ Qo
    if M.dim(a - 2, b + 1):
        Qi = M.Q(a - 2, b + 1).toarray()
        L += Qi @ Qi.conj().T
        wB = np.linalg.eigvalsh(Qi @ Qi.conj().T)
    return np.linalg.eigvalsh((L + L.conj().T) / 2), wB


def test_charge_bookkeeping():
    """Q raises s by 1, V by 1/5; E, K, Z are chiral with delta = 3/5, 4/5, 2/5; particle-hole flips s."""
    N = 7
    for a in range(N + 1):
        for b in range(N + 1):
            assert abs(s_eff(a + 2, b - 1, N) - s_eff(a, b, N) - 1) < 1e-12
            assert abs(s_eff(a + 1, b + 1, N) - s_eff(a, b, N) - 0.2) < 1e-12
            assert abs(s_eff(N - a, N - b, N) + s_eff(a, b, N)) < 1e-12
            assert abs(s_eff(a + 1, b - 1, N) - s_eff(a, b, N) - DELTA["E"]) < 1e-12
            assert abs(s_eff(a + 2, b, N) - s_eff(a, b, N) - DELTA["K"]) < 1e-12
            assert abs(s_eff(a + 1, b, N) - s_eff(a, b, N) - DELTA["Z"]) < 1e-12
    assert abs(lam(10) - 1.0) < 1e-12
    assert bps_fraction(4, 3, 7) > 0 and bps_fraction(0, 0, 7) == 0.0


def test_dressing_targets_follow_the_slot_equations():
    """Targets of the dressing terms (D8 sec. 7) and their chord charges s' + {+-7, +-9, +-11}/10."""
    Np = 8
    expect = {"J10[E h00]": ((1, -1), 0.7), "J10[K h11]": ((1, -1), 0.7), "J10[Z h01]": ((1, -1), 0.7),
              "J01[K+ h00]": ((-2, 0), -0.7), "J01[E+ h11]": ((-2, 0), -0.7), "J01[Z+ h10]": ((-2, 0), -0.7),
              "z11[E h01]": ((1, -2), 0.9), "z00[K h01]": ((2, -1), 1.1),
              "x00[E+ h10]": ((-2, 1), -0.9), "x11[K+ h10]": ((-3, 0), -1.1)}
    for A in range(Np + 1):
        for B in range(Np + 1):
            sn = s_eff(A, B, Np)
            for name, t in dressing_targets(A, B, Np).items():
                (dA, dB), ds = expect[name]
                assert tuple(t["target"]) == (A + dA, B + dB), name
                assert abs(t["s_target"] - (sn + ds)) < 1e-12, name
            res = resonant_terms(A, B, Np)
            if res:
                assert min(abs(abs(sn) - x) for x in (0.2, 0.4, 0.6)) < 1e-12
    assert set(resonant_terms(3, 3, 8)) == {"J10[E h00]", "J10[K h11]", "J10[Z h01]"}     # s' = -0.2
    assert set(resonant_terms(5, 5, 8)) == {"J01[K+ h00]", "J01[E+ h11]", "J01[Z+ h10]"}  # s' = +0.2
    assert resonant_terms(4, 4, 8) == [] and resonant_terms(4, 4, 7) == []


def test_window_contains_every_nontower_bps_sector():
    """N' = 6: every sector with BPS states lies in |s| <= 1/2, except tower sectors (n, n); both models."""
    Np = 6
    for gen in (False, True):
        M = TwoFlavor((random_generic_couplings if gen else random_symmetric_couplings)(Np, np.random.default_rng(3)))
        for a in range(Np + 1):
            for b in range(Np + 1):
                if M.dim(a, b) == 0:
                    continue
                w, _ = laplacian_spectrum(M, a, b)
                nb = int(np.sum(w < 1e-9 * max(1, w.max())))
                if nb and a != b:
                    assert abs(s_eff(a, b, Np)) <= 0.5 + 1e-12, (a, b)


def test_tower_is_the_only_bps_state_outside_the_window():
    """N' = 9, sector (1,1) (s = -0.7): one BPS state (the tower) with symmetrized couplings, none for generic ones."""
    Np = 9
    for gen, expect in ((False, 1), (True, 0)):
        M = TwoFlavor((random_generic_couplings if gen else random_symmetric_couplings)(Np, np.random.default_rng(1)))
        w, _ = laplacian_spectrum(M, 1, 1)
        assert int(np.sum(w < 1e-9 * max(1, w.max()))) == expect
    assert abs(s_eff(1, 1, Np) + 0.7) < 1e-12


def test_window_edge_gap_collapses():
    """Old N = 7: the exact (B) family of the edge sector (4,2) (s = +1/2) reaches down to ~1e-4 of its mean energy,
    while the central sector (4,4) (s = +0.1) has a gap above 1e-2 of its mean."""
    M = TwoFlavor(random_symmetric_couplings(8, np.random.default_rng(1))[:7, :7, :7].copy())
    _, wB = laplacian_spectrum(M, 4, 2)
    wB = wB[wB > 1e-9 * wB.max()]
    assert abs(s_eff(4, 2, 7) - 0.5) < 1e-12 and wB.min() / wB.mean() < 1e-3
    w, wB2 = laplacian_spectrum(M, 4, 4)
    wB2 = wB2[wB2 > 1e-9 * wB2.max()]
    assert wB2.min() / w.mean() > 1e-2


def test_obstruction_weight_matches_chord_point_mass():
    """BPS -> BPS weight of Z from the BPS states of old sector (3,4) at N = 6 (s_out = -0.2, s_in = +0.2): exact
    (exact average over unit-variance couplings) vs the Al-Salam-Chihara point mass (BLY 4.14), within 10%."""
    N = 6
    M = TwoFlavor(random_symmetric_couplings(N + 1, np.random.default_rng(1))[:N, :N, :N].copy())

    def bps(a, b):
        d = M.dim(a, b)
        L = np.zeros((d, d), complex)
        for (x, y), out in (((a, b), True), ((a - 2, b + 1), False)):
            if M.dim(x, y) and M.dim(x + 2, y - 1):
                Qm = M.Q(x, y).toarray()
                L += Qm.conj().T @ Qm if out else Qm @ Qm.conj().T
        w, U = np.linalg.eigh((L + L.conj().T) / 2)
        return U[:, w < 1e-9 * max(1, w.max())]

    W, Wt = bps(3, 4), bps(4, 4)
    num = den = 0.0
    for i in range(N):                                 # exact average: sum over the N single-fermion strings
        X = M.op_block([(1.0, [(M.psi_pos[i], True)])], (3, 4), (4, 4)).toarray() @ W
        num += np.sum(np.abs(Wt.conj().T @ X) ** 2)
        den += np.sum(np.abs(X) ** 2)
    ch = two_point("Z", +1, 3, 4, N)
    assert abs(ch["s_in"] - 0.2) < 1e-12
    assert abs(num / den - ch["point_mass"]) < 0.1 * ch["point_mass"]


if __name__ == "__main__":
    import time
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            t0 = time.time()
            fn()
            print(f"ok  {name}  [{time.time() - t0:.1f}s]", flush=True)
