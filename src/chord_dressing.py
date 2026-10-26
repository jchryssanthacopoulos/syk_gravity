r"""
chord_dressing.py -- the first-order dressing spectrum as a matter-chord problem (docs/derivations.md D7) and the
finite-N quantities used to test it.

First-order dressing operators (progress report 2026-10-13, Prop. 5; src/uplift.py):
    h family:  D  = P E^dag G E P   on the old BPS space B^F,      G = H0^{-1} on exact forms of degree F+2,
    g family:  D' = P E G' E^dag P  on B^{F-1},                    G' = H0^{-1} on coexact forms of degree F-3,
with E = eta ^ = sum_I eta_I psi_I^dag (pairs I of old modes, E|eta_I|^2 = 2).  First-order decoder eigenvalues:
y = ((1 - s)/(1 + s))^2 for s in spec D, spec D'; W^(1) = <y(1-y)>, T2^(1) = <y>.

Chord (double-scaled) prediction.  Wick pairings of the eta's are matter-chord diagrams.  For uncrossed pairings the
chord rules put one energy on every face; faces touching zero-temperature (BPS) segments sit on the unique zero mode of
their charge, so the free cumulants of D are t_k = (1/b) Tr S^k with S = G^{1/2} V G^{1/2}, V = E_eta[E P E^dag],
diagonal in energy with eigenvalue g = (b c0/D_in) A(phi) / (kappa Lambda(phi)).  Hence D is free compound Poisson with
rate r = e/b and jump law g distributed by the B-family density (functions chord_*).

Charge assignments (chord charge s = (degree - N/2)/p): h family s_out = (F - N/2)/p, s_in = s_out + 2/p; g family
mirrored (s -> -s): s_out = -((F-1) - N/2)/p, s_in = s_out + 2/p.  Energy scale kappa fixed by the exact first moment
of H0 in the inner charge sector, 2[C(d,p) + C(N-d,p)] / (q^s + q^-s); c0 = exact tau = 0 value of the two-point
function (2 C(N-F, 2) resp. 2 C(F-1, 2)).  lambda = 2 p^2 / N with N = N' - 1 (old model).  No fitted parameters.
"""
from __future__ import annotations

from math import comb

import numpy as np

from chord_chaos import Sector, WedgeStructure, random_couplings
from dssyk_n2 import b_family_density, bps_continuum_kernel, chord_energy, time_integrated_two_point
from uplift import bps_dim, coexact_forms_dim, exact_forms_dim, split_couplings


# ----------------------------------------------------------------------------------------------------- chord side --

def chord_sectors(Np: int, fam: str, p: int = 3, F: int | None = None) -> dict:
    """Chord parameters of the h ('h') or g ('g') dressing problem at N' (old model N = N' - 1).
    F = N'//2 by default (D7).  Any other F is a window sector (D10 §4).  The resonant families (s_in = 1/2) are the
    h family at l' = 2F - N' = -2 and the g family at l' = +2.  There the jump law has a 1/phi^2 tail."""
    N = Np - 1
    F = Np // 2 if F is None else F
    lam = 2 * p * p / N
    q = np.exp(-lam)
    delta = 2 / p
    if fam == "h":
        s_out, deg_b, deg_in = (F - N / 2) / p, F, F + 2
        c0, b, e = 2 * comb(N - F, 2), bps_dim(N, p, F), exact_forms_dim(N, F + 2, p)
    elif fam == "g":
        s_out, deg_b, deg_in = -((F - 1) - N / 2) / p, F - 1, F - 3
        c0, b, e = 2 * comb(F - 1, 2), bps_dim(N, p, F - 1), coexact_forms_dim(N, F - 3, p)
    else:
        raise ValueError(fam)
    s_in = s_out + delta
    kappa = 2 * (comb(deg_in, p) + comb(N - deg_in, p)) / (q ** s_in + q ** -s_in)
    return dict(Np=Np, N=N, F=F, p=p, lam=lam, q=q, delta=delta, s_out=s_out, s_in=s_in, deg_b=deg_b, deg_in=deg_in,
                c0=c0, kappa=kappa, Din=comb(N, deg_in), b=b, e=e, r=e / b)


def chord_t1(Np: int, fam: str, F: int | None = None) -> float:
    """Mean first-order dressing cost (1/b) Tr D-bar = c0 int_0^inf G(tau) dtau, closed form (2phi1).
    Finite only for s_in > 1/2 (gapped target)."""
    P = chord_sectors(Np, fam, F=F)
    return P["c0"] / P["kappa"] * time_integrated_two_point(P["lam"], P["s_out"], P["delta"])


def chord_jump_law(Np: int, fam: str, nphi: int = 4000, F: int | None = None):
    """Gauss-Legendre nodes in phi: S-eigenvalues g_k and weights w_k (sum 1) of the chord jump law; returns (g, w, P)."""
    P = chord_sectors(Np, fam, F=F)
    x, wg = np.polynomial.legendre.leggauss(nphi)
    phi = (x + 1) * np.pi / 2
    rho = b_family_density(P["lam"], P["s_in"], phi) * wg * np.pi / 2
    E = P["kappa"] * chord_energy(P["lam"], P["s_in"], phi)
    g = P["b"] * P["c0"] / P["Din"] * bps_continuum_kernel(P["lam"], P["s_out"], P["delta"], phi) / E
    return g, rho / rho.sum(), P


def chord_population(Np: int, fam: str, n: int | None = None, nphi: int = 200001, F: int | None = None):
    """n S-eigenvalues (default: the exact count e) placed at the quantiles of the B-family density.
    In a resonant sector (s_in = 1/2, D10 §4) the eigenvalues are unbounded near phi = 0. The largest one is fixed by
    the top quantile (n - 1/2)/n, which plays the role of the finite-n cutoff."""
    P = chord_sectors(Np, fam, F=F)
    n = P["e"] if n is None else n
    phi = np.linspace(0, np.pi, nphi)[1:-1]
    rho = b_family_density(P["lam"], P["s_in"], phi)
    cdf = np.concatenate([[0.0], np.cumsum((rho[1:] + rho[:-1]) / 2 * np.diff(phi))])
    cdf /= cdf[-1]
    ph = np.interp((np.arange(n) + 0.5) / n, cdf, phi)
    E = P["kappa"] * chord_energy(P["lam"], P["s_in"], ph)
    g = P["b"] * P["c0"] / P["Din"] * bps_continuum_kernel(P["lam"], P["s_out"], P["delta"], ph) / E
    return g, P


def wy(s):
    """First-order decoder variables: (w, y) with y = ((1-s)/(1+s))^2, w = y(1-y)."""
    s = np.clip(np.asarray(s, dtype=float), 0.0, None)
    y = ((1 - s) / (1 + s)) ** 2
    return y * (1 - y), y


def separable_samples(sev, cev, b: int, draws: int, rng):
    """Spectra of X^dag X with X = S^1/2 Z C^1/2 (Z: e x b iid complex Gaussian, variance 1/b); sev: S eigenvalues
    (e of them), cev: C eigenvalues (b of them, mean 1; all ones = Wishart / factorized BPS faces)."""
    ss = np.sqrt(np.clip(np.asarray(sev, float), 0, None))
    cc = np.sqrt(np.clip(np.asarray(cev, float), 0, None))
    out = []
    for _ in range(draws):
        Z = (rng.standard_normal((len(ss), b)) + 1j * rng.standard_normal((len(ss), b))) / np.sqrt(2 * b)
        X = (ss[:, None] * Z) * cc[None, :]
        out.append(np.linalg.eigvalsh(X.conj().T @ X))
    return np.concatenate(out)


def chord_prediction(Np: int, rng, b_eff: int = 1500, draws: int = 3, emax: int = 12000) -> dict:
    """b -> infinity chord prediction of W^(1), T2^(1) (b-weighted over the two families), by Monte Carlo of the free
    compound Poisson law at fixed r = e/b with a large effective b (finite-b effects checked to be < 1e-4 in W)."""
    out = {}
    for fam in ("h", "g"):
        g, w, P = chord_jump_law(Np, fam)
        r = P["r"]
        be = b_eff
        if r * be > emax:
            be = max(50, int(emax / r))
        e_eff = int(round(r * be))
        o = np.argsort(g)
        c = np.cumsum(w[o])
        c /= c[-1]
        gp = g[o][np.minimum(np.searchsorted(c, (np.arange(e_eff) + 0.5) / e_eff), len(g) - 1)]
        sp = separable_samples(gp, np.ones(be), be, draws, rng)
        W, y = wy(sp)
        out[fam] = dict(W=float(W.mean()), T2=float(y.mean()), t1=float(r * np.sum(w * g)), t2=float(r * np.sum(w * g * g)),
                        t1_closed=chord_t1(Np, fam), r=r, b=P["b"], b_eff=be)
    bh, bg = out["h"]["b"], out["g"]["b"]
    out["W"] = (bh * out["h"]["W"] + bg * out["g"]["W"]) / (bh + bg)
    out["T2"] = (bh * out["h"]["T2"] + bg * out["g"]["T2"]) / (bh + bg)
    return out


# ------------------------------------------------------------------------------------------------ finite-N side --

def old_model_lean(Np: int, seed: int, p: int = 3):
    """Couplings and old BPS bases exactly as in uplift.OldModel(Np, c, rng) with rng = default_rng(seed) (same couplings
    as results/data/uplift_2026-10-12), without OldModel's dense copies of Q0 and E.  Returns dict."""
    from run_ou_fourarc import bps_basis_fast           # scripts/ on sys.path
    rng = np.random.default_rng(seed)
    F = Np // 2
    c = random_couplings(rng, Sector(Np, p, F=F).M)
    cA, eta = split_couplings(Np, p, c)
    N = Np - 1
    VF, _, _ = bps_basis_fast(Sector(N, p, F=F), cA, bps_dim(N, p, F), rng)
    VF1, _, _ = bps_basis_fast(Sector(N, p, F=F - 1), cA, bps_dim(N, p, F - 1), rng)
    return dict(Np=Np, N=N, F=F, p=p, cA=cA, eta=eta, VF=VF, VF1=VF1, seed=seed)


def old_model_from_uplift(om) -> dict:
    """Lean old-model dict (as old_model_lean) from an uplift.OldModel built for any sector F (D10).  VF and VF1 are the
    admissible sources (ker delta), so DressingFamily works in every window sector."""
    return dict(Np=om.Np, N=om.N, F=om.F, p=om.p, cA=om.cA, eta=om.eta, VF=om.VF, VF1=om.VF1, seed=None)


def pair_maps(N: int, deg: int):
    """For each pair I (combinations order): (rows in L^{deg+2}, cols in L^deg, signs) of psi_I^dag (nonzeros only)."""
    npair = comb(N, 2)
    Ws = WedgeStructure(N, deg, 2)
    out = []
    for t in range(npair):
        e = np.zeros(npair)
        e[t] = 1.0
        M = Ws.matrix(e).tocoo()
        nz = M.data != 0                       # the sparse matrices store explicit zeros for the other pairs
        out.append((M.row[nz], M.col[nz], M.data[nz].real))
    return out


class DressingFamily:
    """Finite-N objects of one family: inner energies/eigenvectors (non-zero modes of H0 on the inner degree),
    the BPS basis, the pair maps, and the matter operator for a given eta."""

    def __init__(self, om: dict, fam: str, inner: str = "eigh"):
        """inner = 'eigh' (default, D7): inner energies from the Gram matrix with a 1e-9 relative cut.
        inner = 'svd': from the singular values of Q0 (cut at 1e-10 sigma_max, as uplift.pinv does). This is needed in
        resonant sectors (D10 §4), where the window-edge multiplet has E_min/E_max down to ~1e-11 at N = 15."""
        N, F, p, cA = om["N"], om["F"], om["p"], om["cA"]
        self.fam, self.N, self.F = fam, N, F
        if fam == "h":
            Qs = WedgeStructure(N, F - 1, p).matrix(cA)                  # L^{F-1} -> L^{F+2}
            A = (Qs @ Qs.conj().T).toarray()                             # Q0 Q0^dag on L^{F+2}
            self.Vb, self.deg_b, self.deg_in = om["VF"], F, F + 2
            self.maps = pair_maps(N, F)                                   # psi^dag: L^F (cols) -> L^{F+2} (rows)
            self._Ws = WedgeStructure(N, F, 2)
        else:
            Qs = WedgeStructure(N, F - 3, p).matrix(cA)                  # L^{F-3} -> L^F
            A = (Qs.conj().T @ Qs).toarray()                             # Q0^dag Q0 on L^{F-3}
            self.Vb, self.deg_b, self.deg_in = om["VF1"], F - 1, F - 3
            self.maps = pair_maps(N, F - 3)                               # psi^dag: L^{F-3} (cols) -> L^{F-1} (rows)
            self._Ws = WedgeStructure(N, F - 3, 2)
        if inner == "svd":
            del A
            Qd = Qs.toarray() if fam == "h" else Qs.conj().T.toarray()       # inner degree = row space of Qd
            U, sv, _ = np.linalg.svd(Qd, full_matrices=False)
            del Qd
            keep = sv > 1e-10 * sv.max()
            self.we, self.Ue = sv[keep] ** 2, np.ascontiguousarray(U[:, keep])
            order = np.argsort(self.we)
            self.we, self.Ue = self.we[order], np.ascontiguousarray(self.Ue[:, order])
        else:
            w, U = np.linalg.eigh(A)
            del A
            keep = w > 1e-9 * w.max()
            self.we, self.Ue = w[keep], np.ascontiguousarray(U[:, keep])
        del U
        self.b = self.Vb.shape[1]
        self.e = len(self.we)

    def _bi_ii(self, r, c):
        return (c, r) if self.fam == "h" else (r, c)                     # (BPS-degree index, inner index)

    def matter(self, eta):
        """Sparse matter operator (inner degree x BPS degree): E for h, E^dag for g."""
        if self.fam == "h":
            return self._Ws.matrix(eta)
        return self._Ws.matrix(eta).conj().T.tocsr()

    def X(self, eta):
        """X = G^1/2 (matter) V_b in the inner energy basis (e x b); first-order spectrum = eig(X^dag X)."""
        EV = self.matter(eta) @ self.Vb
        return (self.Ue.conj().T @ EV) / np.sqrt(self.we)[:, None]

    def V_inner(self):
        """V = E_eta[(matter) P (matter)^dag] on the inner degree (D_in x D_in), P = V_b V_b^dag (scatter, exact)."""
        P = self.Vb @ self.Vb.conj().T
        Din = self.Ue.shape[0]
        V = np.zeros((Din, Din), complex)
        for r, c, sg in self.maps:
            bi, ii = self._bi_ii(r, c)
            V[np.ix_(ii, ii)] += 2 * np.outer(sg, sg) * P[np.ix_(bi, bi)]
        return V

    def lift(self, Op):
        """2 sum_I P (matter_I)^dag Op (matter_I) P as a b x b matrix (Op on the inner degree)."""
        Dout = self.Vb.shape[0]
        M = np.zeros((Dout, Dout), complex)
        for r, c, sg in self.maps:
            bi, ii = self._bi_ii(r, c)
            M[np.ix_(bi, bi)] += np.outer(sg, sg) * Op[np.ix_(ii, ii)]
        out = 2 * (self.Vb.conj().T @ (M @ self.Vb))
        return (out + out.conj().T) / 2

    def G_inner(self):
        return (self.Ue / self.we) @ self.Ue.conj().T


def analyze_family(fam_obj: DressingFamily, eta_own, neta: int, rng, draws: int = 3, third_moment: bool = True,
                   Np: int | None = None, Fsec: int | None = None) -> tuple[dict, dict]:
    """Error decomposition of one family (docs/derivations.md D7 sec. 6):
    own/sampled = exact first order (own eta; neta - 1 fresh etas), kron = separable model with exact S and C = Dbar/t1,
    wishart = exact S and C = 1, chord = chord population (exact e, b) and C = 1, chord_resc = same rescaled to exact t1."""
    F = fam_obj
    b = F.b
    V = F.V_inner()
    Vt = F.Ue.conj().T @ V @ F.Ue
    S = (Vt / np.sqrt(F.we)[:, None]) / np.sqrt(F.we)[None, :]
    S = (S + S.conj().T) / 2
    del Vt
    sev = np.linalg.eigvalsh(S)
    diag_frac = float(np.sum(np.abs(np.diag(S)) ** 2) / np.sum(np.abs(S) ** 2))
    del S
    t1 = float(sev.sum() / b)
    G = F.G_inner()
    Dbar = F.lift(G)
    dbev = np.linalg.eigvalsh(Dbar)
    rec = dict(b=b, e=F.e, t1=t1, t2=float(np.sum(sev ** 2) / b), t3=float(np.sum(sev ** 3) / b), diag_frac_S=diag_frac,
               Dbar_rel_sd=float(dbev.std() / dbev.mean()), Dbar_mean=float(dbev.mean()),
               m2_identity=float(np.mean(dbev ** 2) + np.sum(sev ** 2) / b), Emin=float(F.we.min()))
    if third_moment:
        GVG = G @ V @ G
        D2 = F.lift(GVG)
        del GVG
        rec["m3_nc_exact"] = float(np.trace(Dbar @ Dbar @ Dbar).real / b + 3 * np.trace(Dbar @ D2).real / b + np.sum(sev ** 3) / b)
        rec["m3_nc_factorized"] = float(t1 ** 3 + 3 * t1 * np.sum(sev ** 2) / b + np.sum(sev ** 3) / b)
        del D2
    del V, G
    spec = {}
    samp, m2s, m3s = [], [], []
    for k in range(neta):
        eta = eta_own if k == 0 else rng.standard_normal(len(eta_own)) + 1j * rng.standard_normal(len(eta_own))
        X = F.X(eta)
        ev = np.linalg.eigvalsh(X.conj().T @ X)
        samp.append(ev)
        m2s.append(float(np.mean(ev ** 2)))
        m3s.append(float(np.mean(ev ** 3)))
    spec["own"], spec["sampled"] = samp[0], np.concatenate(samp)
    spec["kron"] = separable_samples(sev, dbev / t1, b, draws, rng)
    spec["wishart"] = separable_samples(sev, np.ones(b), b, draws, rng)
    if Np is not None:
        Fs = fam_obj.F if Fsec is None else Fsec
        Fs = None if Fs == Np // 2 else Fs                                # keep the D7 default path bit-identical
        gpop, P = chord_population(Np, F.fam, n=F.e, F=Fs)
        spec["chord"] = separable_samples(gpop, np.ones(b), b, draws, rng)
        if P["s_in"] > 0.5 + 1e-12:                                       # t1 is finite only for a gapped target
            spec["chord_rescaled"] = separable_samples(gpop * t1 / chord_t1(Np, F.fam, F=Fs), np.ones(b), b, draws, rng)
            rec["t1_chord"] = chord_t1(Np, F.fam, F=Fs)
        else:
            rec["t1_chord"] = float("inf")
    rec["m2_sampled"] = float(np.mean(m2s))
    rec["m2_sampled_se"] = float(np.std(m2s) / np.sqrt(len(m2s)))
    rec["m3_sampled"] = float(np.mean(m3s))
    rec["m3_sampled_se"] = float(np.std(m3s) / np.sqrt(len(m3s)))
    for key, x in spec.items():
        W, y = wy(x)
        rec[key] = dict(W=float(W.mean()), T2=float(y.mean()), m1=float(np.mean(x)),
                        atoms=float(np.mean(x < 1e-9 * max(float(np.max(x)), 1e-300))))
    rec["Dbar_eigs"] = dbev.tolist()
    return rec, spec


def crossing_terms(fam_obj: DressingFamily, nsamp: int, rng, npair: int) -> dict:
    """Unbiased estimators of individual Wick terms of E_eta tr D^3 with three independent Gaussian vectors a, b, g:
    crossing pairing (1/b) Tr[(X_g^dag X_a)(X_b^dag X_g)(X_a^dag X_b)]; identity and full-cycle pairings as checks."""
    cr, idt, cy, tot = [], [], [], []
    for _ in range(nsamp):
        Xa, Xb, Xg = (fam_obj.X(rng.standard_normal(npair) + 1j * rng.standard_normal(npair)) for _k in range(3))
        aa, bb, gg = Xa.conj().T @ Xa, Xb.conj().T @ Xb, Xg.conj().T @ Xg
        ab, bg, ga = Xa.conj().T @ Xb, Xb.conj().T @ Xg, Xg.conj().T @ Xa
        b = fam_obj.b
        cr.append(np.trace(ga @ bg @ ab).real / b)
        idt.append(np.trace(aa @ bb @ gg).real / b)
        cy.append(np.trace(ga @ ab @ bg).real / b)
        tot.append(float(np.mean(np.linalg.eigvalsh(aa) ** 3)))
    ms = lambda v: (float(np.mean(v)), float(np.std(v) / np.sqrt(len(v))))      # noqa: E731
    return dict(crossing=ms(cr), identity=ms(idt), cycle=ms(cy), m3=ms(tot))


def eth_function(fam_obj: DressingFamily):
    """Finite-N 'ETH function': inner energies E_n and V_nn = sum_alpha E_eta |<n| matter |alpha>|^2."""
    V = fam_obj.V_inner()
    Vnn = np.einsum("in,ij,jn->n", fam_obj.Ue.conj(), V, fam_obj.Ue).real
    return fam_obj.we, Vnn


def memory_estimate_gb(Np: int, p: int = 3) -> dict:
    """Rough peak memory (GB) of old_model_lean + analyze_family (h family dominates): BPS-basis construction, then
    P = V_b V_b^dag (D_F^2), the lift scratch (D_F^2), inner Ue/V/G (D_in^2), and the sampled X's (e x b)."""
    N, F = Np - 1, Np // 2
    c16 = 16 / 1e9
    DF, DF1, Din = comb(N, F), comb(N, F - 1), comb(N, F + 2)
    bF, bF1 = bps_dim(N, p, F), bps_dim(N, p, F - 1)
    basis = (comb(N, F + p) * DF + DF * comb(N, F - p) + DF * (comb(N, F + p) + comb(N, F - p)) + 2 * DF * bF) * c16
    held = (DF * bF + DF1 * bF1) * c16
    h_peak = (2 * DF * DF + 4 * Din * Din + 3 * Din * bF + 2 * bF * bF) * c16
    return dict(basis=basis, held=held, analysis=h_peak, total=held + max(basis, h_peak))


def own_site_pair_operator(om: dict, i: int):
    """E_i = (iota_i omega) ^ : L^F -> L^{F+2}, the pair operator built from the old couplings C_{i..}; {Q, c_i} = E_i."""
    from itertools import combinations
    N, F = om["N"], om["F"]
    pairs = list(combinations(range(N), 2))
    pidx = {pq: t for t, pq in enumerate(pairs)}
    coef = np.zeros(len(pairs), complex)
    for k, I in enumerate(combinations(range(N), 3)):            # old-model couplings in combinations order
        if i in I:
            pos = I.index(i)
            coef[pidx[tuple(x for x in I if x != i)]] += (-1) ** pos * om["cA"][k]
    return WedgeStructure(N, F, 2).matrix(coef)


def single_fermion_sum_rule(om: dict, Fm: "DressingFamily"):
    """Exact identity on the BPS space of degree F (h family), D7 sec. 7:
        B_c + sum_i P E_i^dag G E_i P = F * 1,   B_c = sum_i P c_i^dag P_{F-1} c_i P,  E_i = {Q, c_i}.
    (c_i alpha is co-closed for harmonic alpha; its coexact part has squared norm <alpha|E_i^dag G E_i|alpha>.)
    Returns (B_c, sum_i P E_i^dag G E_i P)."""
    N, F = om["N"], om["F"]
    Vb = Fm.Vb
    G = Fm.G_inner()
    W1 = WedgeStructure(N, F - 1, 1)
    Bc = np.zeros((Fm.b, Fm.b), complex)
    SE = np.zeros((Fm.b, Fm.b), complex)
    for i in range(N):
        X = own_site_pair_operator(om, i) @ Vb
        SE += X.conj().T @ (G @ X)
        e = np.zeros(N)
        e[i] = 1.0
        A = om["VF1"].conj().T @ (W1.matrix(e).conj().T @ Vb)
        Bc += A.conj().T @ A
    return Bc, SE
