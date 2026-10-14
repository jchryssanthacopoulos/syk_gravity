#!/usr/bin/env python3
r"""
Heterogeneity of the eta-averaged dressing cost over BPS states (docs/derivations.md D7 sec. 7).

Dbar = E_eta D = 2 sum_I P psi_I G psi_I^dag P (h family; g family analogous).  The chord theory predicts Dbar = t1 * 1
(faces with zero-temperature segments carry the unique zero mode of their charge); at finite N its eigenvalues spread.
Diagnostics per realization and family:
  rsd        relative standard deviation of the eigenvalues of Dbar
  tau        rsd of K(tau) = 2 sum_I P psi_I e^{-tau H0} psi_I^dag P at tau * E_med in TAUS  (K(0) = c0 * 1 exactly)
  bins       Dbar = sum_k Dbar_k over inner-energy quantile bins (lowest 5 %, 5-20 %, 20-50 %, 50-100 % of levels):
             share of t1 and of the variance (traceless Frobenius norm^2) per bin, and the coherence
             C = sum_{k != l} <dD_k, dD_l> / sum_k <dD_k, dD_k>  (0: bins fluctuate independently; > 0: aligned)
  pair model least-squares fit Dbar ~ sum_I g_I B_I, B_I = P psi_I psi_I^dag P = P (1 - n_i)(1 - n_j) P (h) or
             P n_i n_j P (g): explained fraction R2 of the traceless part, for fitted g_I and for the 'natural'
             g_I = Tr(P psi_I G psi_I^dag P) / Tr(B_I) (mean inverse energy of the states created from pair I);
             site model with B_i = P (1 - n_i) P (h) or P n_i P (g).
  free       incoherent (free-probability) estimate of the spread of each bin and of the whole: Dbar_k = sum_n M_n/E_n with
             M_n = 2 sum_I P psi_I |n><n| psi_I^dag P (rank <= number of pairs); if the M_n were independently
             Haar-rotated, the eigenvalue variance would be sum_n [tr M_n^2 - (tr M_n)^2]/E_n^2 (free cumulants add),
             tr = Tr/b.  rsd_free is compared with the observed rsd of the same bin.
  probes     traceless operator correlation corr(A, B) = <dA, dB>/(|dA||dB|) of Dbar with K(tau), with
             B_n = sum_i (P n_i P)^2 (BPS two-point operator of the site occupations, decoder-type) and with the
             single-fermion BPS -> BPS operator B_c = sum_i P c_i^dag P' c_i P (h: P' = BPS at F-1) or
             sum_i P c_i P' c_i^dag P (g: P' = BPS at F).
Usage: .venv/bin/python scripts/analyze_bps_heterogeneity.py --Np 12 --seeds 9000 9001
"""
import argparse
import json
import os
import sys
import time
from math import comb

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from chord_chaos import Sector  # noqa: E402
from chord_dressing import DressingFamily, old_model_lean  # noqa: E402
from run_chord_dressing import OUT, git_state  # noqa: E402

TAUS = (0.1, 0.3, 1.0, 3.0, 10.0)
BINS = (0.0, 0.05, 0.20, 0.50, 1.0)


def rsd(M):
    ev = np.linalg.eigvalsh((M + M.conj().T) / 2)
    return float(ev.std() / ev.mean())


def traceless(M):
    b = M.shape[0]
    return M - np.trace(M) / b * np.eye(b)


def occupations(N, deg):
    """0/1 occupation table n_i(A) for the Fock basis of L^deg (Sector ordering, bit i of the mask)."""
    masks = np.array(Sector(N, 3, F=deg).masks)
    return ((masks[:, None] >> np.arange(N)[None, :]) & 1).astype(float)


def op_corr(A, B):
    dA, dB = traceless(A), traceless(B)
    return float(np.real(np.sum(dA.conj() * dB)) / np.sqrt(np.real(np.sum(np.abs(dA) ** 2)) * np.real(np.sum(np.abs(dB) ** 2))))


def probe_operators(om, Fm):
    """B_n = sum_i (V^dag n_i V)^2 and B_c (single-fermion BPS -> BPS weight) on the family's BPS space."""
    from chord_chaos import WedgeStructure
    N, F = om["N"], om["F"]
    Vb = Fm.Vb
    occ = occupations(N, Fm.deg_b)
    Bn = sum((Vb.conj().T * occ[:, i][None, :]) @ Vb @ ((Vb.conj().T * occ[:, i][None, :]) @ Vb) for i in range(N))
    W1 = WedgeStructure(N, F - 1, 1)                          # e_i ^ : L^{F-1} -> L^F
    Bc = np.zeros((Fm.b, Fm.b), complex)
    for i in range(N):
        e = np.zeros(N)
        e[i] = 1.0
        cdag = W1.matrix(e)                                   # c_i^dag (L^{F-1} -> L^F)
        if Fm.fam == "h":
            Ai = om["VF1"].conj().T @ (cdag.conj().T @ Vb)    # P_{F-1} c_i P_F
        else:
            Ai = om["VF"].conj().T @ (cdag @ Vb)              # P_F c_i^dag P_{F-1}
        Bc += Ai.conj().T @ Ai
    return Bn, Bc


def pair_model(Fm, Dbar, Gdiag_pairs):
    """Least-squares fit of the traceless part of Dbar by sum_I g_I B_I (pair model) and sum_i g_i B_i (site model).
    Uses Tr(B_I B_J) = sum_{A,B} |P_AB|^2 h_I(A) h_J(B) with h diagonal in the Fock basis (no b x b x npair storage)."""
    N = Fm.N
    occ = occupations(N, Fm.deg_b)
    hole = 1 - occ if Fm.fam == "h" else occ                       # h: pair must be empty; g: pair must be filled
    pairs = [(i, j) for i in range(N) for j in range(i + 1, N)]
    Hp = np.stack([hole[:, i] * hole[:, j] for i, j in pairs], axis=1)   # D x npair
    Hs = hole                                                            # D x N
    Vb = Fm.Vb
    b = Vb.shape[1]
    P = Vb @ Vb.conj().T
    Pi = np.abs(P) ** 2
    Pdiag = np.real(np.diag(P))
    del P
    dD = np.real(np.einsum("ab,Ab,Aa->A", Dbar, Vb.conj(), Vb)) if False else np.real(np.sum((Vb @ Dbar) * Vb.conj(), axis=1))
    out = {}
    tD = np.trace(Dbar).real
    nD = np.sum(np.abs(Dbar) ** 2) - tD ** 2 / b                          # ||traceless Dbar||^2
    for name, H in (("pair", Hp), ("site", Hs)):
        Gram = H.T @ Pi @ H                                              # Tr(B_I B_J)
        tr = H.T @ Pdiag                                                 # Tr(B_I)
        Gram_tl = Gram - np.outer(tr, tr) / b                            # <B_I - tr, B_J - tr>
        rhs = H.T @ dD - tr * tD / b                                     # <B_I^tl, Dbar^tl>
        g, *_ = np.linalg.lstsq(Gram_tl, rhs, rcond=1e-10)
        expl = float(g @ rhs)                                            # explained ||.||^2 (projection)
        out[name] = dict(R2=expl / nD, g_rel_sd=float(np.std(g) / abs(np.mean(g))) if name == "pair" else None)
        if name == "pair":
            gn = Gdiag_pairs / np.maximum(tr, 1e-300)                    # natural weights (mean inverse energy)
            resid_n = nD - 2 * gn @ rhs + gn @ Gram_tl @ gn
            # best overall scale for the natural weights
            alpha = (gn @ rhs) / (gn @ Gram_tl @ gn)
            out["pair_natural"] = dict(R2=float(1 - (nD - 2 * alpha * gn @ rhs + alpha ** 2 * gn @ Gram_tl @ gn) / nD),
                                       R2_unscaled=float(1 - resid_n / nD), scale=float(alpha),
                                       g_rel_sd=float(np.std(gn) / np.mean(gn)))
    return out


def level_moments(Fm, chunk=256):
    """tr M_n and tr M_n^2 for every inner level n, via the npair x npair Gram matrices of u_I = V_b^dag psi_I |n>."""
    b = Fm.b
    e = Fm.e
    trM, trM2 = np.zeros(e), np.zeros(e)
    for s0 in range(0, e, chunk):
        sl = slice(s0, min(e, s0 + chunk))
        U = []
        for rr, cc, sg in Fm.maps:
            bi, ii = (cc, rr) if Fm.fam == "h" else (rr, cc)
            U.append(Fm.Vb[bi, :].conj().T @ (sg[:, None] * Fm.Ue[ii, sl]))       # b x L
        U = np.array(U)                                                          # npair x b x L
        Gam = np.einsum("iak,jak->kij", U.conj(), U)                            # L x npair x npair
        trM[sl] = 2 * np.real(np.einsum("kii->k", Gam)) / b
        trM2[sl] = 4 * np.sum(np.abs(Gam) ** 2, axis=(1, 2)) / b
    return trM, trM2


def analyze(Np, seed, fams=("h", "g")):
    om = old_model_lean(Np, seed)
    res = dict(Np=Np, seed=seed)
    for fam in fams:
        Fm = DressingFamily(om, fam)
        we, Ue = Fm.we, Fm.Ue
        b = Fm.b
        G = (Ue / we) @ Ue.conj().T
        Dbar = Fm.lift(G)
        r = dict(b=b, e=Fm.e, rsd=rsd(Dbar), t1=float(np.trace(Dbar).real / b))
        Emed = float(np.median(we))
        r["tau"] = {}
        for tt in TAUS:
            K = Fm.lift((Ue * np.exp(-tt / Emed * we)) @ Ue.conj().T)
            r["tau"][str(tt)] = rsd(K)
        order = np.argsort(we)
        parts = []
        for lo, hi in zip(BINS[:-1], BINS[1:]):
            sel = order[int(lo * len(we)):int(hi * len(we))]
            parts.append(Fm.lift((Ue[:, sel] / we[sel]) @ Ue[:, sel].conj().T))
        dl = [traceless(Pk) for Pk in parts]
        Cmat = np.array([[np.real(np.sum(x.conj() * y)) for y in dl] for x in dl])
        trM, trM2 = level_moments(Fm)
        kap = (trM2 - trM ** 2) / we ** 2
        mean = trM / we
        free_bins = []
        for lo, hi in zip(BINS[:-1], BINS[1:]):
            sel = order[int(lo * len(we)):int(hi * len(we))]
            free_bins.append(float(np.sqrt(kap[sel].sum()) / mean[sel].sum()))
        neff = float(mean.sum() ** 2 / np.sum(mean ** 2))
        r["free"] = dict(rsd_free_bins=free_bins, rsd_free_all=float(np.sqrt(kap.sum()) / mean.sum()), n_eff_levels=neff,
                         npair=len(Fm.maps), rank_ratio=float(np.median(trM ** 2 / np.maximum(trM2, 1e-300)) * 1.0))
        r["bins"] = dict(edges=list(BINS), t1_share=[float(np.trace(Pk).real / np.trace(Dbar).real) for Pk in parts],
                         var_share=[float(Cmat[k, k] / np.trace(Cmat)) for k in range(len(parts))],
                         coherence=float((Cmat.sum() - np.trace(Cmat)) / np.trace(Cmat)),
                         rsd_bin=[rsd(Pk) for Pk in parts], Emax_bin=[float(we[order[int(hi * len(we)) - 1]]) for hi in BINS[1:]])
        # natural pair weights need Tr(P psi_I G psi_I^dag P) per pair
        Gd = []
        for rr, cc, sg in Fm.maps:
            bi, ii = (cc, rr) if fam == "h" else (rr, cc)
            Zb = Fm.Vb[bi, :] * sg[:, None]
            Gd.append(float(np.real(np.sum(Zb.conj() * (G[np.ix_(ii, ii)] @ Zb)))))
        r["models"] = pair_model(Fm, Dbar, np.array(Gd))
        Bn, Bc = probe_operators(om, Fm)
        r["probes"] = dict(Bn=op_corr(Dbar, Bn), Bc=op_corr(Dbar, Bc), Bn_Bc=op_corr(Bn, Bc), rsd_Bn=rsd(Bn), rsd_Bc=rsd(Bc),
                           K={str(tt): op_corr(Dbar, Fm.lift((Ue * np.exp(-tt / Emed * we)) @ Ue.conj().T)) for tt in TAUS})
        res[fam] = r
        del Fm, G, Dbar, parts, dl
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Np", type=int, required=True)
    ap.add_argument("--seeds", type=int, nargs="+", default=[9000])
    ap.add_argument("--fams", nargs="+", default=["h", "g"])
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, f"heterogeneity_N{a.Np}.jsonl")
    with open(os.path.join(a.out, f"meta_heterogeneity_N{a.Np}_{time.strftime('%Y%m%dT%H%M%S')}.json"), "w") as fh:
        json.dump(dict(args=vars(a), git=git_state(), taus=TAUS, bins=BINS), fh, indent=1)
    for seed in a.seeds:
        t0 = time.time()
        res = analyze(a.Np, seed, a.fams)
        res["secs"] = time.time() - t0
        with open(path, "a") as fh:
            fh.write(json.dumps(res) + "\n")
        for fam in a.fams:
            r = res[fam]
            m = r["models"]
            print(f"N'={a.Np} seed={seed} {fam}: rsd {r['rsd']:.3f} | tau: " + " ".join(f"{v:.3f}" for v in r["tau"].values())
                  + f" | bins t1-share {np.round(r['bins']['t1_share'], 3).tolist()} var-share {np.round(r['bins']['var_share'], 3).tolist()}"
                  f" coherence {r['bins']['coherence']:+.3f} | R2 pair {m['pair']['R2']:.3f} natural {m['pair_natural']['R2']:.3f} "
                  f"(scale {m['pair_natural']['scale']:.2f}) site {m['site']['R2']:.3f} [{res['secs']:.0f}s]", flush=True)
            fr = r["free"]
            print(f"      free (incoherent) estimate: bins {np.round(fr['rsd_free_bins'], 3).tolist()} vs observed "
                  f"{np.round(r['bins']['rsd_bin'], 3).tolist()}; whole {fr['rsd_free_all']:.3f} vs {r['rsd']:.3f}; "
                  f"N_eff levels {fr['n_eff_levels']:.1f} of {r['e']}", flush=True)
            pr = r["probes"]
            print(f"      probes: corr(Dbar, B_n) {pr['Bn']:+.3f}  corr(Dbar, B_c) {pr['Bc']:+.3f}  corr(B_n, B_c) {pr['Bn_Bc']:+.3f}  "
                  f"rsd B_n {pr['rsd_Bn']:.3f} B_c {pr['rsd_Bc']:.3f} | corr(Dbar, K(tau)): " + " ".join(f"{v:+.2f}" for v in pr["K"].values()), flush=True)


if __name__ == "__main__":
    main()
