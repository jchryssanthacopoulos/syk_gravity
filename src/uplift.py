r"""
uplift.py -- exact uplift of BPS states of one-flavor N=2 SYK from N = N'-1 to N' modes (docs/derivations.md D6).

Split the added mode e = N'-1 (the last mode) from the old modes o = 0..N'-2:
    Lambda^F(N') = Lambda^F(o) (+) Lambda^{F-1}(o) ^ e ,     Psi = a + b ^ e,
    omega = omega_0 + eta ^ e     (omega_0: couplings not touching e;  eta: (p-1)-form of couplings touching e).
With Q0 = omega_0 ^ and E = eta ^ (E commutes with Q0 for odd p), the BPS conditions on sector F are
    Q0 a = 0,  Q0^dag b = 0,  Q0 b = -eps E a,  Q0^dag a = +eps E^dag b,     eps = (-1)^F,
and their solution is the exact uplift formula
    a = (1 + K'K)^{-1} (h + eps K' g),   b = g - eps K a,        (h, g) in B_o^F (+) B_o^{F-1}  (old harmonic states),
    K = Q0^+ E   (Lambda^F(o) -> coexact Lambda^{F-1}(o)),   K' = (Q0^dag)^+ E^dag   (Lambda^{F-1}(o) -> exact Lambda^F(o)).
First order in eta: decoder eigenvalues 1/(1+sigma^2) (sigma: singular values of K on B_o^F) and sigma'^2/(1+sigma'^2)
(sigma': singular values of K' on B_o^{F-1}).  Exact walls: ker(eta^) on B_o^F, ker(iota_eta) on B_o^{F-1}.
Dressing cost: sigma^2 = eigenvalues of P_B E^dag H0^{-1} E P_B (H0 = old Hamiltonian on Lambda^{F+2}(o)), i.e. a
BPS-projected, Euclidean-time-integrated two-point function of the pair operator E; likewise sigma'^2 with E^dag.
Orthogonal splitting: the uplifts of B_o^F (b coexact) and of B_o^{F-1} (a exact) are exactly orthogonal, so the
decoder is a 2x2 block operator on B = B_h (+) B_g (decoder_blocks).

Conventions follow src/chord_chaos.py (coupling order = itertools.combinations; Fock basis = subset_index masks).
"""
from __future__ import annotations

from math import comb

import numpy as np
import scipy.linalg as sla

from chord_chaos import Sector, WedgeStructure, p_subsets


def bps_dim(N: int, p: int, F: int) -> int:
    """Exact BPS dimension (BBNR eq. 5.25): sum_n (-1)^n C(N, F + n p)."""
    return sum((-1) ** abs(n) * comb(N, F + n * p) for n in range(-(F // p) - 1, (N - F) // p + 2)
               if 0 <= F + n * p <= N)


def exact_forms_dim(N: int, k: int, p: int = 3) -> int:
    """dim of Q0-exact forms in Lambda^k(N), assuming no cohomology in degrees >= k: e_k = C(N,k) - e_{k+p}."""
    return 0 if k > N else comb(N, k) - exact_forms_dim(N, k + p, p)


def coexact_forms_dim(N: int, k: int, p: int = 3) -> int:
    """dim of Q0^dag-exact (coexact) forms in Lambda^k(N), assuming no cohomology in degrees <= k."""
    return 0 if k < 0 else comb(N, k) - coexact_forms_dim(N, k - p, p)


def atom_bound(Np: int, p: int = 3, F: int | None = None):
    """Refined forced-atom count of the decoder (per wall) and the h-wall dimension ratio r = e_{F+2}/b_F."""
    F = Np // 2 if F is None else F
    N = Np - 1
    bh, bg = bps_dim(N, p, F), bps_dim(N, p, F - 1)
    th, tg = exact_forms_dim(N, F + p - 1, p), coexact_forms_dim(N, F - p, p)
    return dict(h_dim=bh, h_target=th, h_atoms=max(0, bh - th), g_dim=bg, g_target=tg, g_atoms=max(0, bg - tg),
                ratio_h=th / bh, ratio_g=tg / bg)


def split_couplings(Np: int, p: int, c: np.ndarray):
    """Old-model couplings (p-subsets not containing the last mode, combinations order of N'-1 modes) and the
    (p-1)-form eta (p-subsets containing the last mode, ordered as combinations of the remaining p-1 old modes)."""
    subs, _ = p_subsets(Np, p)
    last = Np - 1
    cA = np.array([c[k] for k, I in enumerate(subs) if last not in I])
    eta = np.array([c[k] for k, I in enumerate(subs) if last in I])
    return cA, eta


def pinv(A: np.ndarray, rtol: float = 1e-10) -> np.ndarray:
    U, s, Vh = np.linalg.svd(A, full_matrices=False)
    keep = s > rtol * s.max()
    return (Vh[keep].conj().T / s[keep]) @ U[:, keep].conj().T


def pinv_apply(A: np.ndarray, rhs: np.ndarray, rtol: float = 1e-10) -> np.ndarray:
    U, s, Vh = np.linalg.svd(A, full_matrices=False)
    keep = s > rtol * s.max()
    return Vh[keep].conj().T @ ((U[:, keep].conj().T @ rhs) / s[keep][:, None])


class OldModel:
    """Old-model data needed for the uplift at N' = N+1, charge F: harmonic bases B_o^F, B_o^{F-1}, and the maps."""

    def __init__(self, Np: int, c: np.ndarray, rng, p: int = 3, basis_fn=None):
        from run_ou_fourarc import bps_basis_fast  # noqa: WPS433  (scripts/ on sys.path)
        basis_fn = basis_fn or bps_basis_fast
        self.Np, self.p, self.N, self.F = Np, p, Np - 1, Np // 2
        N, F = self.N, self.F
        self.eps = (-1) ** F
        self.cA, self.eta = split_couplings(Np, p, c)
        self.SF, self.SF1 = Sector(N, p, F=F), Sector(N, p, F=F - 1)
        self.VF, _, _ = basis_fn(self.SF, self.cA, bps_dim(N, p, F), rng)
        self.VF1, _, _ = basis_fn(self.SF1, self.cA, bps_dim(N, p, F - 1), rng)
        self.Q0_Fm1 = WedgeStructure(N, F - 1, p).matrix(self.cA).toarray()        # L^{F-1} -> L^{F+2}
        self.Q0_Fm3 = WedgeStructure(N, F - p, p).matrix(self.cA).toarray()        # L^{F-3} -> L^F
        self.E_F = WedgeStructure(N, F, p - 1).matrix(self.eta).toarray()          # L^F -> L^{F+2}
        self.E_Fm3 = WedgeStructure(N, F - p, p - 1).matrix(self.eta).toarray()    # L^{F-3} -> L^{F-1}

    def first_order_singular_values(self):
        """Singular values of K on B_o^F and of K' on B_o^{F-1} (padded with zeros to the full dimensions)."""
        sig = np.linalg.svd(pinv_apply(self.Q0_Fm1, self.E_F @ self.VF), compute_uv=False)
        sigp = np.linalg.svd(pinv_apply(self.Q0_Fm3.conj().T, self.E_Fm3.conj().T @ self.VF1), compute_uv=False)
        sig = np.concatenate([sig, np.zeros(self.VF.shape[1] - len(sig))])
        sigp = np.concatenate([sigp, np.zeros(self.VF1.shape[1] - len(sigp))])
        return sig, sigp

    def first_order_decoder(self):
        sig, sigp = self.first_order_singular_values()
        return np.concatenate([1.0 / (1.0 + sig ** 2), sigp ** 2 / (1.0 + sigp ** 2)])

    def uplift_components(self):
        """Exact uplift (a, b) of the old harmonic bases: columns of (a, b) are Psi(h_j) for h_j in VF, then Psi(g_j) for
        g_j in VF1.  a lives in Lambda^F(o) (n_e = 0), b in Lambda^{F-1}(o) (n_e = 1).  Not orthonormal."""
        K = pinv(self.Q0_Fm1) @ self.E_F                               # L^F -> L^{F-1}
        Kp = pinv(self.Q0_Fm3.conj().T) @ self.E_Fm3.conj().T          # L^{F-1} -> L^F
        M = np.eye(self.SF.D) + Kp @ K
        a = np.linalg.solve(M, np.hstack([self.VF, self.eps * (Kp @ self.VF1)]))
        b = np.hstack([np.zeros((self.SF1.D, self.VF.shape[1])), self.VF1]) - self.eps * (K @ a)
        return a, b

    def decoder_blocks(self):
        """Block form of the decoder mu-operator P U P (U = (-1)^{n_e}) on B = B_h (+) B_g (orthogonal, D6 §5).
        Returns M_h, M_g, C2 in orthonormal bases of the two families (mu = [[M_h, C2], [C2^dag, M_g]], C2 = 2<a_h, a_g>)
        and the max overlap |<Psi_h, Psi_g>| (zero by the splitting theorem)."""
        a, b = self.uplift_components()
        nh = self.VF.shape[1]
        fam = []
        for sl in (slice(0, nh), slice(nh, None)):
            A, B = a[:, sl], b[:, sl]
            w, U = np.linalg.eigh(A.conj().T @ A + B.conj().T @ B)
            S = U / np.sqrt(w)
            fam.append((A @ S, B @ S))
        (Ah, Bh), (Ag, Bg) = fam
        overlap = float(np.abs(Ah.conj().T @ Ag + Bh.conj().T @ Bg).max())
        Mh = Ah.conj().T @ Ah - Bh.conj().T @ Bh
        Mg = Ag.conj().T @ Ag - Bg.conj().T @ Bg
        C2 = Ah.conj().T @ Ag - Bh.conj().T @ Bg
        return Mh, Mg, C2, overlap

    def three_operators(self):
        """K_H = K|_{B_o^F} (into coexact Lambda^{F-1}), K'_H = K'|_{B_o^{F-1}} (into exact Lambda^F), and
        L = P_Co K P_Ex (exact Lambda^F -> coexact Lambda^{F-1}); P_Ex K' P_Co = L^dag (D6 §5).  Returned as dense
        matrices in the old Fock bases: K_H (D_{F-1} x b_F), K'_H (D_F x b_{F-1}), L (D_{F-1} x D_F)."""
        K = pinv(self.Q0_Fm1) @ self.E_F
        Kp = pinv(self.Q0_Fm3.conj().T) @ self.E_Fm3.conj().T
        P_ex = self.Q0_Fm3 @ pinv(self.Q0_Fm3)
        return K @ self.VF, Kp @ self.VF1, K @ P_ex

    def decoder_closed_form(self):
        """Exact decoder as the generalized eigenproblem N v = mu G v on B^0 = B_o^F (+) B_o^{F-1} (D6 §5):
            G = diag(1 + K_H^dag (1+R)^-1 K_H,  1 + K'_H^dag (1+R')^-1 K'_H),
            N = [[1 - K_H^dag (1-R)(1+R)^-2 K_H,  -2 eps K_H^dag L (1+R')^-2 K'_H], [h.c., -1 + K'_H^dag (1-R')(1+R')^-2 K'_H]],
        R = L L^dag, R' = L^dag L.  Returns (N, G); L -> 0 gives the first-order decoder."""
        KH, KpH, L = self.three_operators()
        R, Rp = L @ L.conj().T, L.conj().T @ L
        I1, I0 = np.eye(R.shape[0]), np.eye(Rp.shape[0])
        A1 = np.linalg.solve(I1 + R, KH)                    # (1+R)^-1 K_H
        A0 = np.linalg.solve(I0 + Rp, KpH)                  # (1+R')^-1 K'_H
        nh, ng = KH.shape[1], KpH.shape[1]
        G = np.zeros((nh + ng, nh + ng), complex)
        N = np.zeros_like(G)
        G[:nh, :nh] = np.eye(nh) + KH.conj().T @ A1
        G[nh:, nh:] = np.eye(ng) + KpH.conj().T @ A0
        N[:nh, :nh] = np.eye(nh) - A1.conj().T @ (I1 - R) @ A1
        N[nh:, nh:] = -np.eye(ng) + A0.conj().T @ (I0 - Rp) @ A0
        N[:nh, nh:] = -2 * self.eps * A1.conj().T @ L @ A0
        N[nh:, :nh] = N[:nh, nh:].conj().T
        return N, G

    def exact_uplift(self, SN: Sector):
        """Columns = uplifted BPS states Psi(h_j), Psi(g_j) of the N'-mode model in SN's Fock basis (not orthonormal)."""
        a, b = self.uplift_components()
        idxF = {m: k for k, m in enumerate(self.SF.masks)}
        idxF1 = {m: k for k, m in enumerate(self.SF1.masks)}
        bit = 1 << (self.Np - 1)
        Psi = np.zeros((SN.D, a.shape[1]), complex)
        for k, m in enumerate(SN.masks):
            Psi[k] = b[idxF1[m & ~bit]] if m & bit else a[idxF[m]]
        return Psi


def decoder_from_basis(SN: Sector, V: np.ndarray, Np: int) -> np.ndarray:
    """Decoder eigenvalues lambda (weight in the n_e = 0 sector) for an orthonormal BPS basis V."""
    u = 1.0 - 2.0 * np.array([(m >> (Np - 1)) & 1 for m in SN.masks], float)
    mu = sla.eigvalsh(V.conj().T @ (u[:, None] * V), check_finite=False)
    return (1.0 + mu) / 2.0


def wall_fractions(lam: np.ndarray, eps_list=(0.1, 0.03, 0.01, 0.001)):
    w = np.minimum(lam, 1.0 - lam)
    return [float(np.mean(w < e)) for e in eps_list]
