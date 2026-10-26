"""Tests for the uplift anywhere in the window (src/uplift.py, OldModel(F=...); docs/derivations.md D10).
Plain-script runner: .venv/bin/python tests/test_uplift_window.py"""
import os
import sys

import numpy as np
import scipy.linalg as sla

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from chord_chaos import Sector, random_couplings  # noqa: E402
from uplift import (EXCEPTIONAL_Q3, OldModel, bps_dim_q3, connecting_form, decoder_from_basis, harmonic_basis,  # noqa: E402
                    level, split_couplings, strip_walk_count, uplift_bookkeeping)


def _model(Np, seed, p=3):
    rng = np.random.default_rng(seed)
    return rng, random_couplings(rng, Sector(Np, p, F=Np // 2).M)


def test_bookkeeping_regimes():
    """q = 3: the dressing target of the n_e = 0 slot is l' + 5. It is the window edge at l' = -2 (resonant), inside the
    window at l' = -3 (obstructed), and outside it for l' >= -1 (gapped). The g family mirrors this."""
    reg = {lp: uplift_bookkeeping(12 + (lp % 2), (12 + (lp % 2) + lp) // 2) for lp in range(-3, 4)}
    assert [reg[lp]["h"]["regime"] for lp in range(-3, 4)] == \
        ["obstructed", "resonant", "gapped", "gapped", "gapped", "gapped", "empty"]
    assert [reg[lp]["g"]["regime"] for lp in range(-3, 4)] == \
        ["empty", "gapped", "gapped", "gapped", "gapped", "resonant", "obstructed"]
    b = uplift_bookkeeping(12, 6)
    assert (b["h"]["s_out"], b["h"]["s_in"]) == (1 / 6, 5 / 6)      # D7, even N' half filling
    b5 = uplift_bookkeeping(12, 5, p=5)                             # q = 5: resonant at |l'| = q - 1 = 4
    assert uplift_bookkeeping(12, 4, p=5)["h"]["regime"] == "resonant" and b5["h"]["regime"] == "gapped"


def test_strip_walks_and_bps_dims():
    """CV eq. (130): strip-confined walks = D~_p - D~_{N-p-q}; q = 3 dims 3^{N/2-1}, 2 3^{N/2-1}, 3^{(N-1)/2}."""
    for N in range(4, 15):
        for F in range(N + 1):
            l = level(N, F)
            Dt = lambda P: sum((-1) ** n * __import__("math").comb(N, P - 3 * n) for n in range(P // 3 + 1))  # noqa: E731
            assert strip_walk_count(N, 3, F) == (max(Dt(F) - Dt(N - F - 3), 0) if abs(l) <= 2 else 0)
    assert bps_dim_q3(12, 6) == 486 and bps_dim_q3(12, 5) == 243 and bps_dim_q3(13, 5) == 1 and bps_dim_q3(9, 3) == 3


def test_general_uplift_spans_every_window_sector():
    """For every sector |l'| <= 4 at N' = 5..10: the admissible dimensions sum to the direct BPS dimension, which equals
    the known count (strip walks + exceptional states); the uplift lies in the BPS space and has full rank; the
    closed-form decoder equals the direct one."""
    for Np in range(5, 11):
        rng, c = _model(Np, 9000 + Np)
        for F in range(1, Np):
            if abs(level(Np, F)) > 4:
                continue
            om = OldModel(Np, c, rng, F=F)
            SN = Sector(Np, 3, F=F)
            V = harmonic_basis(SN, c, rng)
            d = om.VF.shape[1] + om.VF1.shape[1]
            assert d == V.shape[1] == bps_dim_q3(Np, F), (Np, F)
            if d:
                Psi = om.exact_uplift(SN)
                sv = np.linalg.svd(Psi, compute_uv=False)
                assert np.linalg.norm(Psi - V @ (V.conj().T @ Psi)) < 1e-12 * np.linalg.norm(Psi)
                assert np.sum(sv > 1e-8 * sv[0]) == d
                Nm, G = om.decoder_closed_form()
                lam = np.sort((1 + sla.eigvalsh(Nm, G)) / 2)
                assert np.abs(lam - np.sort(decoder_from_basis(SN, V, Np))).max() < 1e-10


def test_connecting_map_ranks():
    """Ranks of delta (D10 §3), q = 3:
    - l' = -2 at N = N'-1 = 1 mod 4: delta_h is surjective onto the e exceptional states at l = +3, and the e exceptional
      b-states at l = -3 uplift freely;
    - l' = -3: rank = 3^{N/2-1} - e(N');
    - l' = -4: delta_h is injective on the exceptional states at l = -3;
    - central sectors: no obstruction."""
    for Np in (6, 8, 9, 10, 11):
        rng, c = _model(Np, 7 * Np)
        N = Np - 1
        e = EXCEPTIONAL_Q3.get(N, 0)
        for F in range(1, Np):
            lp = level(Np, F)
            if abs(lp) > 4:
                continue
            ob = OldModel(Np, c, rng, F=F).obstruction
            rh, rg = ob["h"]["rank"], ob["g"]["rank"]
            if lp in (-2, -4):
                assert (rh, rg) == (e, 0), (Np, lp)
            elif lp in (2, 4):
                assert (rh, rg) == (0, e), (Np, lp)
            elif lp == -3:
                assert (rh, rg) == (3 ** (N // 2 - 1) - EXCEPTIONAL_Q3.get(Np, 0), 0), (Np, lp)
            elif lp == 3:
                assert (rh, rg) == (0, 3 ** (N // 2 - 1) - EXCEPTIONAL_Q3.get(Np, 0)), (Np, lp)
            else:
                assert (rh, rg) == (0, 0), (Np, lp)


def test_exceptional_states_are_the_radical_of_the_connecting_form():
    """B(beta, alpha) = int beta ^ eta ^ alpha on H(N, -2) obeys B^T = (-1)^{N/2-1} B. It is antisymmetric for
    N = 0 mod 4 on the odd dimension 3^{N/2-1}, and its radical is the exceptional space at N' = N+1 (1, 3, 1 at N' = 5, 9, 13).
    For N = 2 mod 4 it is symmetric and nondegenerate."""
    for N in (4, 6, 8, 10, 12):
        rng, c = _model(N + 1, 31 + N)
        cA, eta = split_couplings(N + 1, 3, c)
        k = (N - 2) // 2
        M, V, _ = connecting_form(N, 3, k, cA, eta, rng)
        assert M.shape == (3 ** k, 3 ** k)
        assert np.abs(M.T - (-1) ** k * M).max() < 1e-12 * np.abs(M).max()
        s = np.linalg.svd(M, compute_uv=False)
        assert np.sum(s < 1e-9 * s[0]) == EXCEPTIONAL_Q3.get(N + 1, 0), N


def test_mirror_symmetry_and_splitting_off_center():
    """Particle-hole: spec lambda(F) = 1 - spec lambda(N' - F) per realization; the orthogonal splitting B = B_h + B_g
    holds off half filling (N' = 10, l' = -2 has both families: 78 ordinary h-states, 3 exceptional g-states)."""
    for Np in (10, 11):
        rng, c = _model(Np, 5)
        for F in range(1, Np // 2 + 1):
            if abs(level(Np, F)) > 3:
                continue
            lams = []
            for FF in (F, Np - F):
                om = OldModel(Np, c, rng, F=FF)
                if om.VF.shape[1] + om.VF1.shape[1] == 0:
                    lams.append(np.zeros(0))
                    continue
                Mh, Mg, C2, overlap = om.decoder_blocks()
                assert overlap < 1e-12
                Nm, G = om.decoder_closed_form()
                lams.append(np.sort((1 + sla.eigvalsh(Nm, G)) / 2))
            assert len(lams[0]) == len(lams[1]), (Np, F)
            if len(lams[0]):
                assert np.abs(lams[0] - np.sort(1 - lams[1])).max() < 1e-10, (Np, F)
    om = OldModel(10, _model(10, 5)[1], np.random.default_rng(0), F=4)
    assert (om.VF.shape[1], om.VF1.shape[1]) == (78, 3)


def test_dyson_form_of_the_self_dressed_metric():
    """D10 §4: K_H^dag (1+R)^-1 K_H = e^dag (S^2 + W)^-1 e with W = U^dag E P_Ex E^dag U (exact). W = S R S in the singular
    basis, so W_jj = s_j^2 <v_j|R|v_j>, and the self-energy vanishes on the near-zero inner modes."""
    rng, c = _model(12, 9000)
    om = OldModel(12, c, rng, F=5)                       # resonant h family (target = window edge l = +3)
    Nm, G = om.decoder_closed_form()
    dy = om.dyson("h")
    d = dy["e"].shape[1]
    Gd = np.eye(d) + dy["e"].conj().T @ np.linalg.solve(np.diag(dy["s2"]) + dy["W"], dy["e"])
    assert np.abs(Gd - G[:d, :d]).max() < 1e-10 * np.abs(G[:d, :d]).max()
    r = np.real(np.diag(dy["W"])) / dy["s2"]
    assert r.max() < 1.0 and dy["s2"].min() / dy["s2"].mean() < 1e-4   # hard edge of the edge multiplet; r bounded


def test_kanazawa_wettig_appendix_b_rows():
    """KW 2017 App. B (q = 3): N+_f = rank(Q: L^f -> L^{f+3}) (coexact), N-_f = rank(Q: L^{f-3} -> L^f) (exact),
    Nz_f = C(N, f) - N+_f - N-_f (BPS), for one random realization.  Rows copied from the KW tables."""
    from math import comb

    from chord_chaos import WedgeStructure
    kw = {5: ([1, 4, 1, 0, 0, 0], [0, 1, 9, 9, 1, 0]),
          9: ([1, 9, 36, 80, 36, 9, 1, 0, 0, 0], [0, 0, 0, 3, 81, 81, 3, 0, 0, 0]),
          13: ([1, 13, 78, 285, 702, 1208, 702, 285, 78, 13, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 729, 729, 1, 0, 0, 0, 0, 0])}
    for N, (Nplus, Nz) in kw.items():
        c = random_couplings(np.random.default_rng(N), comb(N, 3))

        def rank(k):
            W = WedgeStructure(N, k, 3)
            if not (W.shape[0] and W.shape[1]):
                return 0
            s = np.linalg.svd(W.matrix(c).toarray(), compute_uv=False)
            return int(np.sum(s > 1e-9 * s[0]))
        plus = [rank(f) for f in range(N + 1)]
        minus = [plus[f - 3] if f >= 3 else 0 for f in range(N + 1)]
        assert plus == Nplus and minus == [0, 0, 0] + Nplus[:-3], N
        assert [comb(N, f) - plus[f] - minus[f] for f in range(N + 1)] == Nz, N


def test_resonant_chord_weight_and_jump_tail():
    """D10 §4 (chord limit, s_in = 1/2): the ASC weight is finite at phi = 0, f(0) = 2 (q^2;q^2) / (pi (b;q^2)),
    b = q^{1/2 + delta - s_out}.  The jump law g ~ G0/phi^2, so g phi^2 tends to a constant: the K^{-1/2} tail."""
    from chord_dressing import chord_sectors
    from dssyk_n2 import bps_continuum_kernel, bps_continuum_weight, chord_energy, poch
    for Np in (12, 16):
        P = chord_sectors(Np, "h", F=(Np - 2) // 2)
        lam, s_out, delta = P["lam"], P["s_out"], P["delta"]
        assert abs(P["s_in"] - 0.5) < 1e-12
        q = np.exp(-lam)
        f0 = 2 * poch(q * q, q * q) / (np.pi * poch(q ** (0.5 + delta - s_out), q * q))
        assert abs(float(bps_continuum_weight(lam, s_out, delta, np.array([1e-4]))[0]) - f0) < 1e-6 * f0
        phi = np.array([1e-2, 1e-3, 1e-4])
        g_phi2 = bps_continuum_kernel(lam, s_out, delta, phi) / chord_energy(lam, P["s_in"], phi) * phi ** 2
        assert np.abs(g_phi2 / g_phi2[-1] - 1).max() < 1e-3


if __name__ == "__main__":
    import time
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            t0 = time.time()
            fn()
            print(f"ok  {name}  [{time.time() - t0:.1f}s]", flush=True)
