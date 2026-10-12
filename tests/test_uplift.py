"""Tests for src/uplift.py (docs/derivations.md D6)."""
import os
import sys

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))

from chord_chaos import Sector, random_couplings  # noqa: E402
from run_ou_fourarc import bps_basis_fast  # noqa: E402
from uplift import OldModel, atom_bound, bps_dim, decoder_from_basis  # noqa: E402


def _exact_basis(Np, c, rng):
    SN = Sector(Np, 3, F=Np // 2)
    V, _, _ = bps_basis_fast(SN, c, bps_dim(Np, 3, Np // 2), rng)
    return SN, V


def test_exact_uplift_formula_spans_bps_space():
    for Np in (8, 9, 10):
        rng = np.random.default_rng(Np)
        SN = Sector(Np, 3, F=Np // 2)
        c = random_couplings(rng, SN.M)
        SN, V = _exact_basis(Np, c, rng)
        Psi = OldModel(Np, c, rng).exact_uplift(SN)
        res = np.linalg.norm(Psi - V @ (V.conj().T @ Psi)) / np.linalg.norm(Psi)
        sv = np.linalg.svd(Psi, compute_uv=False)
        assert res < 1e-12 and (sv > 1e-8 * sv[0]).sum() == V.shape[1]


def test_first_order_exact_atoms_match_decoder_atoms():
    # exact walls (lambda in {0,1}) are kernels of eta^ / iota_eta on old BPS states: counted identically by
    # first order (sigma = 0) and by the exact decoder
    Np = 11
    rng = np.random.default_rng(3)
    SN = Sector(Np, 3, F=Np // 2)
    c = random_couplings(rng, SN.M)
    lam1 = OldModel(Np, c, rng).first_order_decoder()
    SN, V = _exact_basis(Np, c, rng)
    lam = decoder_from_basis(SN, V, Np)
    n1 = int(np.sum(np.minimum(lam1, 1 - lam1) < 1e-9))
    ne = int(np.sum(np.minimum(lam, 1 - lam) < 1e-9))
    assert n1 == ne and ne >= atom_bound(Np)['h_atoms'] + atom_bound(Np)['g_atoms']


def test_orthogonal_splitting_and_block_decoder():
    # D6 §5: <Psi(h), Psi(g)> = 0 exactly; the block operator [[M_h, C2], [C2^dag, M_g]] reproduces the exact decoder
    for Np in (9, 10, 11):
        rng = np.random.default_rng(100 + Np)
        SN = Sector(Np, 3, F=Np // 2)
        c = random_couplings(rng, SN.M)
        om = OldModel(Np, c, rng)
        Mh, Mg, C2, overlap = om.decoder_blocks()
        assert overlap < 1e-12
        mu_blocks = np.linalg.eigvalsh(np.block([[Mh, C2], [C2.conj().T, Mg]]))
        SN, V = _exact_basis(Np, c, rng)
        mu_exact = 2 * decoder_from_basis(SN, V, Np) - 1
        assert np.abs(np.sort(mu_blocks) - np.sort(mu_exact)).max() < 1e-10


def test_three_operator_closed_form():
    # D6 §5: P_Ex K' P_Co = L^dag (anti-Hermitian coupling), and the closed-form (N, G) pencil gives the exact decoder;
    # self-dressing is monotone (ordered h-eigenvalues of the diagonal block >= first order)
    import scipy.linalg as sla
    from uplift import pinv
    for Np in (10, 11):
        rng = np.random.default_rng(200 + Np)
        SN = Sector(Np, 3, F=Np // 2)
        c = random_couplings(rng, SN.M)
        om = OldModel(Np, c, rng)
        KH, KpH, L = om.three_operators()
        Kp = pinv(om.Q0_Fm3.conj().T) @ om.E_Fm3.conj().T
        P_ex = om.Q0_Fm3 @ pinv(om.Q0_Fm3)
        P_co = pinv(om.Q0_Fm1) @ om.Q0_Fm1
        assert np.abs(P_ex @ Kp @ P_co - L.conj().T).max() < 1e-12 * max(1.0, np.abs(L).max())
        N, G = om.decoder_closed_form()
        mu = sla.eigvalsh(N, G)
        SN, V = _exact_basis(Np, c, rng)
        assert np.abs(np.sort(mu) - np.sort(2 * decoder_from_basis(SN, V, Np) - 1)).max() < 1e-10
        nh = KH.shape[1]
        mu_self = np.sort(sla.eigvalsh(N[:nh, :nh], G[:nh, :nh]))
        sig, _ = om.first_order_singular_values()
        assert np.all(mu_self >= np.sort((1 - sig ** 2) / (1 + sig ** 2)) - 1e-10)


def test_dressing_cost_is_inverse_hamiltonian_two_point_function():
    # D6 §3: sigma^2 = spec of P_B E^dag H0^{-1} E P_B, H0 = Q0 Q0^dag + Q0^dag Q0 on Lambda^{F+2}(o)
    from chord_chaos import WedgeStructure
    Np = 12                                                             # N = 11: F + 5 = 11 <= N, both terms
    rng = np.random.default_rng(7)
    c = random_couplings(rng, Sector(Np, 3, F=Np // 2).M)
    om = OldModel(Np, c, rng)
    sig, _ = om.first_order_singular_values()
    N, F = om.N, om.F
    Qin = om.Q0_Fm1                                                     # L^{F-1} -> L^{F+2}
    H0 = Qin @ Qin.conj().T
    if F + 5 <= N:
        Qout = WedgeStructure(N, F + 2, 3).matrix(om.cA).toarray()      # L^{F+2} -> L^{F+5}
        H0 = H0 + Qout.conj().T @ Qout
    EV = om.E_F @ om.VF
    Dh = EV.conj().T @ np.linalg.pinv(H0, hermitian=True) @ EV
    assert np.allclose(np.sort(np.linalg.eigvalsh(Dh)), np.sort(sig ** 2), atol=1e-9)


def test_atom_bound_values():
    assert atom_bound(10)['h_atoms'] == 45 and atom_bound(12)['h_atoms'] == 79
    assert atom_bound(14)['h_atoms'] == 27 and atom_bound(15)['h_atoms'] == 0 and atom_bound(16)['h_atoms'] == 0


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
