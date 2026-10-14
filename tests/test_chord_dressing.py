"""Tests for the BPS -> continuum chord closed forms (src/dssyk_n2.py, D7), the transfer-matrix engine
(src/chord_transfer.py) and the finite-N dressing machinery (src/chord_dressing.py).
Run: .venv/bin/python tests/test_chord_dressing.py"""
import os
import sys
from math import comb

import numpy as np
from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
from chord_transfer import ChordN2, bps_matter_spectral_weights, geig, ground_state_amplitudes, qpoch  # noqa: E402
from dssyk_n2 import (bps_bps_weight, bps_continuum_kernel, bps_continuum_weight, bps_fraction_chord,  # noqa: E402
                      b_family_density, chord_energy, log_poch, time_integrated_two_point)


def qhermite(x, Q, M):
    H = np.zeros(M + 1)
    H[0] = 1.0
    if M >= 1:
        H[1] = 2 * x
    for m in range(1, M):
        H[m + 1] = 2 * x * H[m] - (1 - Q ** m) * H[m - 1]
    return H


def test_ground_state_is_q_hermite_and_engine_zero_mode():
    lam, s = 1.4, 1 / 6
    q = np.exp(-lam)
    x, o = ground_state_amplitudes(q, s, 12)
    H = qhermite(-np.cosh(lam * (s + 0.5)), q * q, 12)
    pred = np.array([q ** (1.5 * m) * H[m] / qpoch(q * q, q * q, m) for m in range(13)])
    assert np.max(np.abs(pred - x)) < 1e-14
    C = ChordN2(lam, s, 60)
    T, G, _ = C.bosonic()
    w, V = geig(T, G)
    assert abs(w[0]) < 1e-12
    Psi = V[:, 0] / V[0, 0]
    states = [C.states[i] for i in C.boson]
    for k, st in enumerate(states[1:9], start=1):
        n = len(st) // 2
        ref = x[n] if st[0] == "X" else o[n]
        assert abs(Psi[k] - ref) < 1e-12, (st, Psi[k], ref)


def test_bbnr_eigenvector_sign():
    """With alpha_n = q^{n/2} H_n(cos phi|q^2)/(q^2;q^2)_n (BBNR 5.10), v = sum alpha_n Q|n+1/2,X> has eigenvalue
    q^-s + q^{s-1} + 2 q^{-1/2} cos phi, i.e. Lambda_s(pi - phi) in the convention of BBNR (5.3)."""
    lam, s, phi = 1.2, 0.7, 0.9
    q = np.exp(-lam)
    C = ChordN2(lam, s, 40)
    H = qhermite(np.cos(phi), q * q, 40)
    v = np.zeros(len(C.states))
    for n in range(19):
        st = tuple("X" if j % 2 == 0 else "O" for j in range(2 * n + 1))
        e = np.zeros(len(C.states))
        e[C.idx[st]] = 1.0
        v += q ** (n / 2) * H[n] / qpoch(q * q, q * q, n) * (C.Q @ e)
    Tv = C.T_full @ v
    lam_expected = chord_energy(lam, s, np.pi - phi)
    sl = slice(0, 20)
    assert np.max(np.abs(Tv[sl] - lam_expected * v[sl])) < 1e-9 * np.max(np.abs(v[sl]))


def test_closed_form_matches_engine_moments():
    for lam, s_out, d in ((1.4, 1 / 6, 2 / 3), (1.8, 0.0, 2 / 3), (0.7, 1 / 6, 2 / 3), (1.1, 0.0, 0.8)):
        wi, f = bps_matter_spectral_weights(lam, s_out, d, lmax=200)
        s_in = s_out + d
        for j in (0, 1, 2, 3, -1, -2):
            a = quad(lambda x: bps_continuum_weight(lam, s_out, d, x) * chord_energy(lam, s_in, x) ** j, 0, np.pi,
                     limit=500, epsabs=1e-14, epsrel=1e-13)[0]
            bb = float(np.sum(f * wi ** float(j)))
            assert abs(a - bb) < 1e-8 * abs(bb), (lam, s_out, d, j, a, bb)


def test_weight_is_density_times_kernel_and_normalized():
    lam, s_out, d = 1.3, 0.1, 0.6
    x = np.linspace(0.05, 3.1, 7)
    assert np.allclose(bps_continuum_weight(lam, s_out, d, x),
                       b_family_density(lam, s_out + d, x) * bps_continuum_kernel(lam, s_out, d, x), rtol=1e-12)
    assert abs(quad(lambda t: bps_continuum_weight(lam, s_out, d, t), 0, np.pi, limit=400)[0] - 1) < 1e-10


def test_point_mass_is_bly_4_14():
    """s_in < 1/2: point mass (q^{1-2 s_in};q^2)/(q^{1-2 s_out};q^2) = 1 - continuum = engine zero-mode weight = BLY (4.14)
    normalized amplitude x sqrt(D(s_in)/D(s_out))."""
    for lam, s_out, d in ((1.0, -0.3, 0.5), (1.5, -0.4, 2 / 3), (1.2, 0.1, 0.2)):
        q = np.exp(-lam)
        s_in = s_out + d
        pm = bps_bps_weight(lam, s_out, d)
        cont = quad(lambda x: bps_continuum_weight(lam, s_out, d, x), 0, np.pi, limit=400)[0]
        assert abs(pm + cont - 1) < 1e-10
        wi, f = bps_matter_spectral_weights(lam, s_out, d, lmax=200)
        assert abs(f[np.argmin(np.abs(wi))] - pm) < 1e-9
        # BLY (4.14) with j_R = -s_out, j'_R = j_R - d, without the q^{2 Delta j_R} prefactor
        jR, jp = -s_out, -s_out - d
        L = lambda a: log_poch(a, q * q)                                                     # noqa: E731
        R = np.exp(0.5 * (L(q ** (1 - 2 * jR)) + L(q ** (1 + 2 * jp)) - L(q ** (1 + 2 * jR)) - L(q ** (1 - 2 * jp))))
        assert abs(R * np.sqrt(bps_fraction_chord(lam, s_in) / bps_fraction_chord(lam, s_out)) - pm) < 1e-10


def test_time_integral_2phi1():
    for lam, s_out, d in ((1.4, 1 / 6, 2 / 3), (0.4, 0.0, 2 / 3), (1.2, -0.2, 0.9)):
        num = quad(lambda x: bps_continuum_weight(lam, s_out, d, x) / chord_energy(lam, s_out + d, x), 0, np.pi,
                   limit=1000, epsabs=1e-14, epsrel=1e-13)[0]
        assert abs(num - time_integrated_two_point(lam, s_out, d)) < 1e-11 * num


def test_finite_N_identities_and_agreement_with_oldmodel():
    """N' = 10: (i) Tr G V = Tr Dbar (both = b t1); (ii) the own-eta spectrum equals OldModel's first-order sigma^2;
    (iii) scatter construction of V equals the dense sum."""
    from chord_dressing import DressingFamily, old_model_lean
    from chord_chaos import Sector, random_couplings
    from uplift import OldModel
    Np, seed = 10, 9000
    om = old_model_lean(Np, seed)
    rng = np.random.default_rng(seed)
    c = random_couplings(rng, Sector(Np, 3, F=Np // 2).M)
    ref = OldModel(Np, c, rng)
    sig, sigp = ref.first_order_singular_values()
    for fam, s_ref in (("h", sig ** 2), ("g", sigp ** 2)):
        Fm = DressingFamily(om, fam)
        X = Fm.X(om["eta"])
        ev = np.sort(np.clip(np.linalg.eigvalsh(X.conj().T @ X), 0, None))
        assert np.max(np.abs(ev - np.sort(s_ref))) < 1e-9 * max(1.0, ev.max())
        V = Fm.V_inner()
        G = Fm.G_inner()
        assert abs(np.trace(G @ V).real - np.trace(Fm.lift(G)).real) < 1e-9 * abs(np.trace(G @ V))
        Vd = np.zeros_like(V)
        npair = comb(om["N"], 2)
        P = Fm.Vb @ Fm.Vb.conj().T
        for t in range(npair):
            e = np.zeros(npair)
            e[t] = 1.0
            M = Fm.matter(e).toarray()
            Vd += 2 * M @ P @ M.conj().T
        assert np.max(np.abs(Vd - V)) < 1e-12


def test_single_fermion_dressing_sum_rule():
    """B_c + sum_i P E_i^dag G E_i P = F * 1 on the BPS space (h family), E_i = {Q, c_i} = (iota_i omega)^."""
    from chord_dressing import DressingFamily, old_model_lean, single_fermion_sum_rule
    for Np in (10, 11):
        om = old_model_lean(Np, 9001)
        Fm = DressingFamily(om, "h")
        Bc, SE = single_fermion_sum_rule(om, Fm)
        assert np.abs(Bc + SE - om["F"] * np.eye(Fm.b)).max() < 1e-11


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
