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

Any sector of the window (docs/derivations.md D10).  OldModel(Np, c, rng, F=F) works for every charge F.
- New-sector level l' = 2F - N'. The sources sit at old levels l'+1 (n_e = 0 slot, degree F) and l'-1 (n_e = 1 slot,
  degree F-1).
- E maps the n_e = 0 slot to level l' + 2p - 1, and E^dag maps the n_e = 1 slot to level l' - 2p + 1.
- When such a target level lies in the window, E h can have a harmonic part. The connecting map
  delta = P_H E P_H (CV eq. 25) then obstructs some sources.
- The formula above holds verbatim on the admissible sources A_h = ker(P_H E|B_o^F) and A_g = ker(P_H E^dag|B_o^{F-1}).
  These are stored as VF, VF1. The full harmonic bases are BF, BF1.
- At half filling the targets lie outside the window, so A = B and nothing changes.
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
    if A.size == 0:
        return np.zeros((A.shape[1], A.shape[0]), dtype=A.dtype)
    U, s, Vh = np.linalg.svd(A, full_matrices=False)
    keep = s > rtol * s.max()
    return (Vh[keep].conj().T / s[keep]) @ U[:, keep].conj().T


def pinv_apply(A: np.ndarray, rhs: np.ndarray, rtol: float = 1e-10) -> np.ndarray:
    if A.size == 0:
        return np.zeros((A.shape[1], rhs.shape[1]), dtype=complex)
    U, s, Vh = np.linalg.svd(A, full_matrices=False)
    keep = s > rtol * s.max()
    return Vh[keep].conj().T @ ((U[:, keep].conj().T @ rhs) / s[keep][:, None])


def range_basis(A: np.ndarray, rtol: float = 1e-10):
    """Orthonormal basis of the column space of A and the nonzero singular values (descending), via a thin SVD."""
    if A.size == 0:
        return np.zeros((A.shape[0], 0), complex), np.zeros(0)
    U, s, _ = np.linalg.svd(A, full_matrices=False)
    keep = s > rtol * s.max()
    return U[:, keep], s[keep]


# --------------------------------------------------------------------------- any sector of the window (D10) --

def level(N: int, F: int) -> int:
    """Level l = 2F - N of degree F in the N-mode model (CV eq. 29)."""
    return 2 * F - N


def strip_walk_count(N: int, p: int, F: int) -> int:
    """Number of +-1 walks of length N from 0 to l = 2F - N that stay in |x| <= p - 1 (CV App. A, S_p).
    Unconditionally a walk count.  It equals the ordinary (non-exceptional) BPS dimension conditional on CV's rank
    formula (Conjecture A.1): S_p = D~_p - D~_{N-p-q} (CV eq. 130)."""
    w = {0: 1}
    for _ in range(N):
        nw = {}
        for x, m in w.items():
            for y in (x - 1, x + 1):
                if abs(y) <= p - 1:
                    nw[y] = nw.get(y, 0) + m
        w = nw
    return w.get(2 * F - N, 0)


EXCEPTIONAL_Q3 = {5: 1, 9: 3, 13: 1, 17: 1}      # per window edge, q = 3, N <= 17 (KW App. B; CV Fig. 1)


def bps_dim_q3(N: int, F: int) -> int:
    """BPS dimension of the q = 3 model, N <= 17: strip walks plus the known exceptional states at l = +-3."""
    if abs(level(N, F)) == 3:
        if N > 17:
            raise ValueError("exceptional counts are known only for N <= 17")
        return EXCEPTIONAL_Q3.get(N, 0)
    return strip_walk_count(N, 3, F)


def uplift_bookkeeping(Np: int, F: int, p: int = 3) -> dict:
    """Levels, degrees and chord charges of the uplift into sector (N', F) (D10 §1).
    Slots: a (n_e = 0, old degree F, level l'+1) and b (n_e = 1, old degree F-1, level l'-1).
    Dressing targets: E a at degree F+p-1 (level l'+2p-1); E^dag b at degree F-p (level l'-2p+1).
    The dressing returns to the other slot via Q0^+ (K a in the b slot, K' b in the a slot).
    Chord charges s = l/(2p): h family (s_out, s_in) = (s_a, s_a + 2/p), g family mirrored (-s_b, -s_b + 2/p).
    Regime of each family, decided by s_in against 1/2:
      'gapped' (s_in > 1/2: the target has no BPS states, its continuum is gapped);
      'resonant' (s_in = 1/2: the target is the window edge |l| = p, where the continuum is gapless);
      'obstructed' (s_in < 1/2: the target holds an extensive BPS space, so delta can be nonzero);
      'empty' (the slot has no BPS states: |level| > p)."""
    N = Np - 1
    lp = level(Np, F)
    la, lb = lp + 1, lp - 1
    ta, tb = la + 2 * (p - 1), lb - 2 * (p - 1)
    out = dict(Np=Np, N=N, F=F, p=p, lp=lp, slot_a=dict(level=la, degree=F), slot_b=dict(level=lb, degree=F - 1),
               target_a=dict(level=ta, degree=F + p - 1), target_b=dict(level=tb, degree=F - p))
    for fam, ls, lt in (("h", la, ta), ("g", -lb, -tb)):
        s_out, s_in = ls / (2 * p), lt / (2 * p)
        if abs(ls) > p:
            reg = "empty"
        elif abs(s_in - 0.5) < 1e-12:
            reg = "resonant"
        elif s_in < 0.5:
            reg = "obstructed"
        else:
            reg = "gapped"
        out[fam] = dict(s_out=s_out, s_in=s_in, regime=reg)
    return out


def harmonic_basis(S: Sector, c: np.ndarray, rng, d_expected: int | None = None, rtol: float = 1e-10) -> np.ndarray:
    """Orthonormal basis (D x d) of the harmonic forms ker Q n ker Q^dag in Lambda^F, for any degree F.
    This includes the window edges, degrees with F < p or F + p > N (where one of the two maps is empty), and degrees
    with no harmonic forms (d = 0).
    - The ranges of Q^dag and Q are taken from thin SVDs, keeping singular values sigma > rtol sigma_max. This is the
      same criterion as pinv and range_basis.
    - The random draws are those of run_ou_fourarc.bps_basis_fast, so the spanned space is the same whenever that
      function's Gram threshold (sigma^2 > 1e-10 sigma_max^2) decides correctly.
    - The Gram threshold fails at the window edge once the edge multiplet's smallest mode drops below 1e-5 sigma_max.
      That happened at N = 15 (D10 §4), so the singular values themselves are used here."""
    if S.D == 0:
        return np.zeros((0, 0), complex)
    has_up, has_down = S.up.shape[0] > 0, S.down.shape[1] > 0
    Rs = [np.zeros((S.D, 0), complex)]
    if has_up:
        Rs.append(range_basis(S.up.matrix(c).toarray().conj().T, rtol)[0])
    if has_down:
        Rs.append(range_basis(S.down.matrix(c).toarray(), rtol)[0])
    R = np.hstack(Rs)
    d = S.D - R.shape[1]
    if d_expected is not None and d != d_expected:
        raise RuntimeError(f"BPS dimension {d} != expected {d_expected}")
    V = rng.standard_normal((S.D, d)) + 1j * rng.standard_normal((S.D, d))
    for _ in range(2):
        V = V - R @ (R.conj().T @ V)
        V, _ = np.linalg.qr(V)
    return V


def wedge_pairing(N: int, k: int) -> np.ndarray:
    """Top-form pairing W (C(N,k) x C(N,N-k)): int beta ^ gamma = beta^T W gamma for beta in Lambda^k,
    gamma in Lambda^{N-k} (Fock bases ordered as chord_chaos.Sector / q_scan.subset_index)."""
    from q_scan import merge_sign, subset_index
    rows, _ = subset_index(N, k)
    _, cmap = subset_index(N, N - k)
    full = (1 << N) - 1
    W = np.zeros((len(rows), len(cmap)))
    for i, S in enumerate(rows):
        W[i, cmap[full ^ S]] = merge_sign(S, full ^ S)
    return W


def connecting_form(N: int, p: int, k: int, cA: np.ndarray, eta: np.ndarray, rng, k2: int | None = None):
    """Bilinear connecting form B(beta, alpha) = int beta ^ eta ^ alpha on harmonic forms of the N-mode model.
    alpha is in H^k and beta in H^{k2}, with k2 = N - k - (p-1) by default. Returns (M, V_k, V_k2) with
    M[i, j] = B(beta_i, alpha_j).
    For k2 = k, i.e. k = (N-p+1)/2 (level 1-p), the form obeys B(alpha, beta) = (-1)^k B(beta, alpha), so
    M^T = (-1)^k M.
    The kernel of M (as a map on alpha) is ker(P_H E P_H | H^k), i.e. the kernel of CV's delta[alpha] = [eta ^ alpha]."""
    k2 = N - k - (p - 1) if k2 is None else k2
    Vk = harmonic_basis(Sector(N, p, F=k), cA, rng)
    Vk2 = Vk if k2 == k else harmonic_basis(Sector(N, p, F=k2), cA, rng)
    E = WedgeStructure(N, k, p - 1).matrix(eta).toarray()            # Lambda^k -> Lambda^{k+p-1} = Lambda^{N-k2}
    M = Vk2.T @ wedge_pairing(N, k2) @ (E @ Vk)
    return M, Vk, Vk2


class OldModel:
    """Old-model data needed for the uplift at N' = N+1, charge F: harmonic bases B_o^F, B_o^{F-1}, and the maps.

    F = None means half filling (F = N'//2), with the original behaviour: same random draws, and the BPS dimensions
    checked against bps_dim.  For any other F (D10), the harmonic bases are computed without an expected dimension.
    The obstruction maps are then
        delta_h = P_H E   : B_o^F     -> H^{F+p-1}   (computed as (1 - P_exact) E h, since E h is Q0-closed),
        delta_g = P_H E^+ : B_o^{F-1} -> H^{F-p}     (computed as (1 - P_coexact) E^dag g, since E^dag g is coclosed).
    The admissible sources are VF = B_o^F n ker delta_h and VF1 = B_o^{F-1} n ker delta_g.  All the methods below
    act on them.  The full harmonic bases are kept as BF, BF1, and the singular values of delta in self.obstruction."""

    def __init__(self, Np: int, c: np.ndarray, rng, p: int = 3, basis_fn=None, F: int | None = None,
                 obstruction_rtol: float = 1e-8):
        from run_ou_fourarc import bps_basis_fast  # noqa: WPS433  (scripts/ on sys.path)
        basis_fn = basis_fn or bps_basis_fast
        half = F is None
        self.Np, self.p, self.N, self.F = Np, p, Np - 1, (Np // 2 if half else F)
        N, F = self.N, self.F
        self.eps = (-1) ** F
        self.cA, self.eta = split_couplings(Np, p, c)
        self.SF, self.SF1 = Sector(N, p, F=F), Sector(N, p, F=F - 1)
        if half:
            self.VF, _, _ = basis_fn(self.SF, self.cA, bps_dim(N, p, F), rng)
            self.VF1, _, _ = basis_fn(self.SF1, self.cA, bps_dim(N, p, F - 1), rng)
        else:
            self.VF = harmonic_basis(self.SF, self.cA, rng)
            self.VF1 = harmonic_basis(self.SF1, self.cA, rng)
        self.BF, self.BF1 = self.VF, self.VF1
        self.Q0_Fm1 = WedgeStructure(N, F - 1, p).matrix(self.cA).toarray()        # L^{F-1} -> L^{F+2}
        self.Q0_Fm3 = WedgeStructure(N, F - p, p).matrix(self.cA).toarray()        # L^{F-3} -> L^F
        self.E_F = WedgeStructure(N, F, p - 1).matrix(self.eta).toarray()          # L^F -> L^{F+2}
        self.E_Fm3 = WedgeStructure(N, F - p, p - 1).matrix(self.eta).toarray()    # L^{F-3} -> L^{F-1}
        self.obstruction = None
        if not half:
            self.obstruction = self._admissible(obstruction_rtol)

    def _admissible(self, rtol: float) -> dict:
        """Restrict VF, VF1 to the kernels of the connecting maps (D10 §2).  Returns the singular values of
        delta_h, delta_g in units of the norm scale of E h (resp. E^dag g), and the ranks."""
        out = {}
        for fam, B, X, Q, side in (("h", self.BF, self.E_F, self.Q0_Fm1, "range"),
                                   ("g", self.BF1, self.E_Fm3.conj().T, self.Q0_Fm3, "corange")):
            if B.shape[1] == 0 or X.shape[0] == 0:
                out[fam] = dict(sv=[], rank=0, dim=int(B.shape[1]), admissible=int(B.shape[1]))
                continue
            Y = X @ B                                                        # E h (closed) / E^dag g (coclosed)
            R, _ = range_basis(Q if side == "range" else Q.conj().T)        # exact (resp. coexact) forms of the target
            Z = Y - R @ (R.conj().T @ Y)                                     # = P_H E h  (resp. P_H E^dag g)
            scale = np.linalg.norm(Y) / np.sqrt(Y.shape[1])
            _, s, Vh = np.linalg.svd(Z, full_matrices=True)
            s = s / scale if scale > 0 else np.zeros_like(s)
            rank = int(np.sum(s > rtol))
            if rank:
                B = B @ Vh[rank:].conj().T
                if fam == "h":
                    self.VF = B
                else:
                    self.VF1 = B
            out[fam] = dict(sv=[float(x) for x in s], rank=rank, dim=int(Y.shape[1]), admissible=int(B.shape[1]))
        return out

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
        Ov = Ah.conj().T @ Ag + Bh.conj().T @ Bg
        overlap = float(np.abs(Ov).max()) if Ov.size else 0.0
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

    def dyson(self, fam: str = "h", rtol: float = 1e-10) -> dict:
        """First-order dressing and its exact resummation in the singular basis of the old supercharge (D10 §4).
        h family:
        - Q0 = sum_j s_j |u_j><v_j| (coexact Lambda^{F-1} -> exact Lambda^{F+p-1}), so K_H h = sum_j v_j <u_j|E h>/s_j.
        - With e = U^dag E V_F, first order is K_H^dag K_H = e^dag S^-2 e.
        - The exact self-dressed metric is K_H^dag (1+R)^-1 K_H = e^dag (S^2 + W)^-1 e, where W = U^dag E P_Ex E^dag U
          (P_Ex: exact forms of the n_e = 0 slot).  This is the propagator of H0 + E P_Ex E^dag on the target exact forms
          in place of H0^-1: a Dyson resummation with self-energy E P_Ex E^dag >= 0.
        g family: the mirror, with Q0^dag (exact Lambda^F <- coexact Lambda^{F-p}), E^dag, and P_Co on Lambda^{F-1}.
        Returns dict(s2 = inner energies s_j^2 (descending), e, W)."""
        if fam == "h":
            Q, X, Vs = self.Q0_Fm1, self.E_F, self.VF
            P_slot, _ = range_basis(self.Q0_Fm3)                    # exact forms of Lambda^F (n_e = 0 slot)
        else:
            Q, X, Vs = self.Q0_Fm3.conj().T, self.E_Fm3.conj().T, self.VF1
            P_slot, _ = range_basis(self.Q0_Fm1.conj().T)           # coexact forms of Lambda^{F-1} (n_e = 1 slot)
        if Q.size == 0:
            return dict(s2=np.zeros(0), e=np.zeros((0, Vs.shape[1]), complex), W=np.zeros((0, 0), complex))
        U, s, _ = np.linalg.svd(Q, full_matrices=False)
        keep = s > rtol * s.max()
        U, s = U[:, keep], s[keep]
        e = U.conj().T @ (X @ Vs)
        Y = U.conj().T @ (X @ P_slot)
        return dict(s2=s ** 2, e=e, W=Y @ Y.conj().T)


def decoder_from_basis(SN: Sector, V: np.ndarray, Np: int) -> np.ndarray:
    """Decoder eigenvalues lambda (weight in the n_e = 0 sector) for an orthonormal BPS basis V."""
    u = 1.0 - 2.0 * np.array([(m >> (Np - 1)) & 1 for m in SN.masks], float)
    mu = sla.eigvalsh(V.conj().T @ (u[:, None] * V), check_finite=False)
    return (1.0 + mu) / 2.0


def wall_fractions(lam: np.ndarray, eps_list=(0.1, 0.03, 0.01, 0.001)):
    w = np.minimum(lam, 1.0 - lam)
    return [float(np.mean(w < e)) for e in eps_list]
