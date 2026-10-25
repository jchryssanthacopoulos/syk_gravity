#!/usr/bin/env python3
r"""
Analysis and figures for the chord picture of the two-flavor uplift (docs/derivations.md D9; data from
scripts/run_two_flavor_chords.py in results/data/two_flavor_chords_2026-10-25).

Checks (printed and written to summary.json in the data directory):
  window     BPS sectors vs the chord window |s| < 1/2 (tower and charge-space-boundary sectors listed); exact BPS
             fraction vs D(s); the gap collapse at the window edge s = +-1/2; N' = 9, 10 boundary sectors
  blocks     chord dictionary: BPS -> BPS (obstruction) weights and time-integrated BPS -> continuum functions, exact vs
             chord, by model, N and distance of the target from the window edge
  resonance  exact dressing kappa = 1/mu - 1 per new sector vs the chord charge s' and the resonance flags; tower
  tower      the tower's individual emissions vs black-hole averages; cancellation residuals
  pairing    slot-weight factorization for fortuitous states vs the tower's pairing

Figures (results/figures/two_flavor_chords_2026-10-25/):
  fig_window       sector maps (BPS fraction, window band, tower diagonal, gapless edges) and fraction vs D(s)
  fig_dictionary   obstruction weights, continuum blocks, and the edge gap vs s, exact vs chord
  fig_resonance    dressing per sector vs s' for N' = 4..8 (both models), resonance positions, tower
  fig_tower        tower vs black holes: individual emissions and pairing correlator

Usage: .venv/bin/python scripts/analyze_two_flavor_chords.py [--no-figs]
"""
import argparse
import glob
import json
import os
import sys
from collections import defaultdict
from math import comb

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from two_flavor_chords import bps_fraction, edge_energies, lam  # noqa: E402

DATA = os.path.join(ROOT, "results/data/two_flavor_chords_2026-10-25")
FIGS = os.path.join(ROOT, "results/figures/two_flavor_chords_2026-10-25")
MODELS = ("sym", "generic")
MLAB = {"sym": "symmetrized", "generic": "generic"}


def load(task, model):
    out = []
    for f in sorted(glob.glob(os.path.join(DATA, f"{task}_*_{model}.jsonl"))):
        out += [json.loads(line) for line in open(f)]
    return out


def gm(x):
    x = np.asarray([v for v in x if np.isfinite(v)])
    return (float(np.exp(x.mean())), float(x.std()), len(x)) if len(x) else (float("nan"), float("nan"), 0)


# ---------------------------------------------------------------------------------------------------- checks
def check_window(model):
    R = load("window", model)
    out = dict(by_Np={})
    for Np in sorted({r["Np"] for r in R}):
        rs = [r for r in R if r["Np"] == Np]
        full = all(r["Np"] <= 8 for r in rs)
        bps = [r for r in rs if r["nbps"] > 0]
        outside = sorted({(r["a"], r["b"], round(r["s"], 2), r["nbps"]) for r in bps if not r["in_window"] and not r["edge"]})
        ratio = [r["frac"] / r["D_chord"] for r in rs if r["in_window"] and r["dim"] >= 20 and r["D_chord"] > 0]
        cen = [r["frac"] / r["D_chord"] for r in rs if r["in_window"] and r["dim"] >= 20 and r["D_chord"] > 0
               and 0.3 * Np <= r["a"] <= 0.7 * Np and 0.3 * Np <= r["b"] <= 0.7 * Np]
        edge = []
        bulk = []
        for r in rs:
            if r["dim"] < 100 or r["m1"] <= 0:
                continue
            if r["edge"]:
                g = r["gapB"] if r["s"] > 0 else r["gapA"]
                if g:
                    edge.append(g / r["m1"])
            else:
                gs = [g for g in (r["gapA"], r["gapB"]) if g]
                if gs:
                    bulk.append(min(gs) / r["m1"])
        out["by_Np"][Np] = dict(full=full, n_sectors=len(rs), n_bps_sectors=len(bps), outside_window=outside,
                                frac_over_D=[float(np.mean(ratio)), float(np.min(ratio)), float(np.max(ratio))] if ratio else None,
                                frac_over_D_central=[float(np.mean(cen)), float(np.min(cen)), float(np.max(cen))] if cen else None,
                                edge_gap_rel=[float(np.min(edge)), float(np.max(edge))] if edge else None,
                                bulk_gap_rel_median=float(np.median(bulk)) if bulk else None,
                                bulk_gap_rel_min=float(np.min(bulk)) if bulk else None)
        if Np >= 9:
            out["by_Np"][Np]["boundary"] = sorted({(r["seed"], r["a"], r["b"], round(r["s"], 2), r["nbps"]) for r in rs})
    return out, R


def check_blocks(model):
    R = load("blocks", model)
    out = {}
    for N in sorted({r["N"] for r in R}):
        rs = [r for r in R if r["N"] == N and r["c0"] > 0 and r["d"] >= 5 and r["T_exact"] > 0]
        pm = [np.log(r["pm_chord"] / r["pm_exact"]) for r in rs if r["pm_chord"] > 0 and r["pm_exact"] and r["pm_exact"] > 0]
        near = [np.log(r["T_chord"] / r["T_exact"]) for r in rs if r["T_chord"] and abs(abs(r["s_in"]) - 0.5) < 0.15]
        far = [np.log(r["T_chord"] / r["T_exact"]) for r in rs if r["T_chord"] and abs(abs(r["s_in"]) - 0.5) >= 0.15]
        tdim = lambda r: comb(N, r["target"][0]) * comb(N, r["target"][1])          # noqa: E731
        edge = [r["T_exact"] / r["c0"] for r in rs if r["edge_target"] and tdim(r) >= 50]   # large (gapless) edge sectors
        other = [r["T_exact"] / r["c0"] for r in rs if not r["edge_target"]]
        out[N] = dict(lam=lam(N), point_mass=gm(pm), cont_near=gm(near), cont_far=gm(far),
                      edge_T_over_c0_median=float(np.median(edge)) if edge else None,
                      other_T_over_c0_median=float(np.median(other)) if other else None, n=len(rs))
    return out, R


def check_resonance(model):
    R = load("resonance", model)
    out = {}
    for N in sorted({r["N"] for r in R}):
        Np = N + 1
        rs = [r for r in R if r["N"] == N and r["d_new"] >= 30]
        Jl = [r for r in rs if resonant_target_dim(r, "J") >= 50]                 # J-resonant, large (gapless) edge target
        Js = [r for r in rs if 0 < resonant_target_dim(r, "J") < 50]              # J-resonant, small boundary edge target
        Z = [r for r in rs if r["resonant"] and resonant_target_dim(r, "J") == 0]
        non = [r for r in rs if not r["resonant"]]

        def stats(lst):
            if not lst:
                return None
            return dict(n=len(lst), sectors=sorted({(r["A"], r["B"]) for r in lst}),
                        mean_kappa_median=float(np.median([r["kappa_mean"] for r in lst])),
                        q50_kappa_median=float(np.median([r["kappa_q"][2] for r in lst])),
                        q90_kappa_median=float(np.median([r["kappa_q"][4] for r in lst])),
                        max_kappa_median=float(np.median([r["kappa_max"] for r in lst])),
                        max_kappa_max=float(np.max([r["kappa_max"] for r in lst])))
        tw = [r.get("tower_bare_fraction") for r in rs if "tower_bare_fraction" in r]
        byD = defaultdict(list)
        for r in rs:
            if not (0 < resonant_target_dim(r, "J") < 50):
                byD[round(current_edge_distance(r["s_new"]), 2)].append(r)
        out[Np] = dict(J_resonant_large_target=stats(Jl), J_resonant_small_target=stats(Js), z_resonant=stats(Z),
                       nonresonant=stats(non),
                       by_current_edge_distance={D: stats(v) for D, v in sorted(byD.items())},
                       tower_bare_fraction_min=float(min(tw)) if tw else None, seeds=sorted({r["seed"] for r in rs}))
    return out, R


def resonant_target_dim(r, kind):
    """Largest dimension of an old-model target sector of a resonant dressing term of kind 'J' (currents) or 'z'
    (single-flavor slots); 0 if there is none.  The chord continuum and its gapless edge need a large target sector:
    small boundary sectors have no near-zero modes (D9 sec. 5)."""
    N = r["N"]
    dims = [comb(N, t["target"][0]) * comb(N, t["target"][1]) for k, t in r["terms"].items()
            if t["resonant"] and (k.startswith("J") if kind == "J" else not k.startswith("J"))]
    return max(dims) if dims else 0


def check_tower():
    R = load("tower", "sym")
    out = dict(n=len(R), resid_max=float(max(r["resid"] for r in R)) if R else None, edge=[], other=[])
    for r in R:
        if r["tower_G_first"] is None or r["bh_G_first"] is None:
            continue
        if abs(r["tower_G_first"] - r["bh_G_first"]) <= 1e-9 * abs(r["bh_G_first"]):
            continue                    # B(n,n) is spanned by the tower alone (n = N): no black-hole comparison
        row = dict(N=r["N"], seed=r["seed"], n=r["n"], current=r["current"], s_target=r["s_target"], gap_rel=r["gap_rel"],
                   tower=r["tower_G_first"], bh=r["bh_G_first"])
        (out["edge"] if abs(abs(r["s_target"]) - 0.5) < 1e-9 else out["other"]).append(row)
    return out, R


def check_twisted():
    """Coherence of black-hole E- and K-emissions (tower projected out), symmetrized vs generic, per (N, n)."""
    S, G = load("twisted", "sym"), load("twisted", "generic")
    out = {}
    for N in sorted({r["N"] for r in S}):
        for n in sorted({r["n"] for r in S if r["N"] == N}):
            a = [r["c2_mean"] for r in S if r["N"] == N and r["n"] == n]
            b = [r["c2_mean"] for r in G if r["N"] == N and r["n"] == n]
            t = [r["tower_cos2"] for r in S if r["N"] == N and r["n"] == n]
            s_t = [r["s_target"] for r in S if r["N"] == N and r["n"] == n][0]
            out[f"N{N}_n{n}"] = dict(s_target=s_t, sym=[float(np.mean(a)), float(np.std(a))],
                                     generic=[float(np.mean(b)), float(np.std(b))] if b else None,
                                     tower_cos2_min=float(min(t)), n_seeds=len(a))
    return out, S, G


def check_pairing(model):
    R = load("pairing", model)
    out = {}
    for N in sorted({r["N"] for r in R}):
        rs = [r for r in R if r["N"] == N and r.get("d_fort", 0) >= 10]
        cs = [r["c_mean"] for r in rs]
        dpsi = [abs(r["n_psi"] - r["n_psi_site"]) for r in rs]
        tw = [(r["A"], r["tower_c"]) for r in R if r["N"] == N and "tower_c" in r]
        out[N + 1] = dict(n=len(rs), c_mean_range=[float(min(cs)), float(max(cs))] if cs else None,
                          c_mean_median=float(np.median(cs)) if cs else None,
                          n_psi_dev_max=float(max(dpsi)) if dpsi else None,
                          tower_c=sorted(set((int(a), round(c, 6)) for a, c in tw)))
    return out, R


# ---------------------------------------------------------------------------------------------------- figures
def fig_window(plt, W):
    fig, axs = plt.subplots(1, 3, figsize=(12.5, 4.0))
    for ax, Np in zip(axs[:2], (7, 8)):
        rs = [r for r in W["sym"] if r["Np"] == Np and r["seed"] == 1]
        M = np.full((Np + 1, Np + 1), np.nan)
        for r in rs:
            M[r["b"], r["a"]] = r["frac"]
        im = ax.imshow(M, origin="lower", cmap="viridis", vmin=0, vmax=1, extent=(-0.5, Np + 0.5, -0.5, Np + 0.5))
        x = np.linspace(-0.5, Np + 0.5, 50)
        for c, ls in ((0, "-"), (2.5, "--"), (-2.5, "--")):
            ax.plot(x, 2 * x - Np / 2 - c, "w" + ls, lw=1.2)                 # 2a - b - N'/2 = 5 s
        ax.plot(x, x, color="tab:red", lw=1, ls=":")
        for r in rs:
            g = r["gapB"] if r["s"] > 0 else r["gapA"]
            if r["edge"] and g and r["m1"] > 0 and g / r["m1"] < 1e-2:
                ax.plot(r["a"], r["b"], marker="s", mfc="none", mec="tab:orange", ms=11, mew=1.8)
            if r["nbps"] == 0:
                ax.text(r["a"], r["b"], "0", ha="center", va="center", fontsize=6, color="0.6")
        ax.set_xlim(-0.5, Np + 0.5)
        ax.set_ylim(-0.5, Np + 0.5)
        ax.set_xlabel(r"$N_\psi$")
        ax.set_ylabel(r"$N_\chi$")
        ax.set_title(f"$N'={Np}$ (symmetrized, seed 1)")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03, label="BPS fraction $d/\\dim$")
    ax = axs[2]
    cmap = plt.get_cmap("plasma")
    for Np in range(5, 9):
        for model, mk in (("sym", "o"), ("generic", "^")):
            rs = [r for r in W[model] if r["Np"] == Np and r["in_window"] and r["dim"] >= 20]
            ax.scatter([r["D_chord"] for r in rs], [r["frac"] for r in rs], s=14, marker=mk, color=cmap((Np - 5) / 4.0),
                       alpha=0.75, label=f"$N'={Np}$" if model == "sym" else None)
    ax.plot([0, 1], [0, 1], "k--", lw=0.8)
    ax.set_xlim(0, 0.8)
    ax.set_ylim(0, 0.8)
    ax.set_xlabel(r"chord BPS fraction $D(s)$ (BLY 3.17)")
    ax.set_ylabel("exact BPS fraction")
    ax.set_title("window sectors, dim $\\geq 20$")
    ax.legend(loc="upper left")
    fig.tight_layout()
    return fig


def fig_dictionary(plt, B, W):
    fig, axs = plt.subplots(1, 3, figsize=(12.5, 3.9))
    ax = axs[0]
    for model, mk in (("sym", "o"), ("generic", "^")):
        rs = [r for r in B[model] if r["pm_chord"] > 0 and r["pm_exact"] and r["d"] >= 5]
        ax.scatter([r["pm_chord"] for r in rs], [r["pm_exact"] for r in rs], s=10, marker=mk, alpha=0.6, label=MLAB[model])
    ax.plot([0, 1], [0, 1], "k--", lw=0.8)
    ax.set_xlabel("chord point mass (BLY 4.14)")
    ax.set_ylabel(r"exact $\|P_H O h\|^2/\|O h\|^2$")
    ax.set_title("obstruction weights, $N=3$-$7$")
    ax.legend()
    ax = axs[1]
    for model, mk in (("sym", "o"), ("generic", "^")):
        rs = [r for r in B[model] if r["T_chord"] and r["T_exact"] > 0 and r["d"] >= 5]
        near = np.array([abs(abs(r["s_in"]) - 0.5) < 0.15 for r in rs])
        x = np.array([r["T_chord"] / r["c0"] for r in rs])
        y = np.array([r["T_exact"] / r["c0"] for r in rs])
        ax.scatter(x[~near], y[~near], s=10, marker=mk, color="tab:blue", alpha=0.6,
                   label=r"$|s_{in}\mp\frac{1}{2}|\geq0.15$" if model == "sym" else None)
        ax.scatter(x[near], y[near], s=10, marker=mk, color="tab:orange", alpha=0.7,
                   label=r"$|s_{in}\mp\frac{1}{2}|=0.1$" if model == "sym" else None)
        ed = [r for r in B[model] if r["edge_target"] and r["d"] >= 5 and r["T_exact"] > 0]
        ax.scatter([3.0] * len(ed), [r["T_exact"] / r["c0"] for r in ed], s=12, marker=mk, color="tab:red", alpha=0.7,
                   label=r"$s_{in}=\pm\frac{1}{2}$ (chord $=\infty$)" if model == "sym" else None)
    lo, hi = 1e-3, 3
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, 4)
    ax.set_xlabel(r"chord $\frac{1}{\kappa}\int f/\Lambda$")
    ax.set_ylabel(r"exact $\langle O^\dagger G O\rangle_{BPS}/c_0$")
    ax.set_title("time-integrated BPS$\\to$continuum")
    ax.legend(fontsize=7, loc="upper left")
    ax = axs[2]
    for Np, col in ((6, "tab:green"), (7, "tab:purple"), (8, "tab:brown")):
        rs = [r for r in W["sym"] if r["Np"] == Np and r["seed"] == 1 and r["dim"] >= 50 and r["m1"] > 0]
        for r in rs:
            for fam, g, mk in (("B", r["gapB"], "o"), ("A", r["gapA"], "s")):
                if g:
                    ax.scatter(r["s"], g / r["m1"], s=12, marker=mk, color=col, alpha=0.7)
        ss = np.linspace(-1.2, 1.2, 400)
        L = lam(Np)
        qv = np.exp(-L)
        eB = 2 * qv ** -0.5 * (np.cosh(L * (ss - 0.5)) - 1) / (qv ** ss + qv ** -ss)
        eA = 2 * qv ** -0.5 * (np.cosh(L * (ss + 0.5)) - 1) / (qv ** ss + qv ** -ss)
        ax.plot(ss, np.maximum(eB, 1e-7), color=col, lw=1, label=f"$N'={Np}$")
        ax.plot(ss, np.maximum(eA, 1e-7), color=col, lw=1, ls="--")
    ax.set_yscale("log")
    ax.set_ylim(1e-5, 3)
    ax.axvline(0.5, color="0.5", lw=0.6)
    ax.axvline(-0.5, color="0.5", lw=0.6)
    ax.set_xlabel("chord charge $s$")
    ax.set_ylabel(r"lowest nonzero eigenvalue / mean")
    ax.set_title("continuum edges: exact vs chord")
    ax.legend(fontsize=7)
    fig.tight_layout()
    return fig


def current_edge_distance(s_new):
    """Distance of the uplift currents' targets (s' + 7/10 for J10, raising; s' - 7/10 for J01, lowering) from their
    gapless edges (+1/2, -1/2): min(|s' + 1/5|, |s' - 1/5|)."""
    return min(abs(s_new + 0.2), abs(s_new - 0.2))


def fig_resonance(plt, Rz):
    Nps = [6, 7, 8]
    fig = plt.figure(figsize=(14, 3.8))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 1.25])
    axs = [fig.add_subplot(gs[0, i]) for i in range(3)]
    for ax, Np in zip(axs, Nps):
        for model, off, col in (("sym", -0.02, "tab:blue"), ("generic", 0.02, "tab:green")):
            rs = [r for r in Rz[model] if r["N"] == Np - 1 and r["d_new"] >= 4]
            bysec = defaultdict(list)
            for r in rs:
                bysec[(r["A"], r["B"], round(r["s_new"], 2))].append(r)
            for (A, B, s), lst in bysec.items():
                k = np.concatenate([np.array(r["kappa"]) for r in lst])
                k = k[k > 1e-10]
                if not len(k):
                    continue
                x = s + off + 0.012 * (B - A) / Np
                ax.plot([x, x], [np.quantile(k, 0.5), np.quantile(k, 0.9)], color=col, lw=1, alpha=0.6)
                ax.scatter(x, np.mean(k), marker="x", color=col, s=14)
                ax.scatter(x, np.quantile(k, 0.5), marker="o", color=col, s=8)
                if A == B and model == "sym":
                    ax.scatter(x, 1.3e-3, marker="*", color="tab:red", s=40, zorder=5)
        for sres in (-0.4, -0.2, 0.2, 0.4):
            ax.axvline(sres, color="tab:red" if abs(sres) == 0.2 else "tab:orange", lw=0.8, ls=":")
        ax.set_yscale("log")
        ax.set_ylim(1e-3, 3e3)
        ax.set_xlim(-0.65, 0.65)
        ax.set_xlabel(r"chord charge $s'$ of the new sector")
        ax.set_title(f"$N'={Np}$")
    axs[0].set_ylabel(r"dressing $\kappa=1/\mu-1$")
    axs[0].scatter([], [], marker="x", color="k", label="mean")
    axs[0].scatter([], [], marker="o", color="k", s=8, label="median (bar: to 90%)")
    axs[0].scatter([], [], marker="*", color="tab:red", label=r"tower: $\kappa=0$")
    axs[0].plot([], [], color="tab:blue", label="symmetrized")
    axs[0].plot([], [], color="tab:green", label="generic")
    axs[0].legend(loc="upper left", fontsize=6.5)
    ax = fig.add_subplot(gs[0, 3])
    for model, mk, col in (("sym", "o", "tab:blue"), ("generic", "^", "tab:green")):
        small = [r for r in Rz[model] if r["d_new"] >= 30 and r["N"] >= 4 and 0 < resonant_target_dim(r, "J") < 50]
        ax.scatter([current_edge_distance(r["s_new"]) + (0.006 if model == "generic" else -0.006) for r in small],
                   [r["kappa_max"] for r in small], s=12, marker=mk, color="0.6", alpha=0.7,
                   label="resonant, small boundary target" if model == "sym" else None)
        rs = [r for r in Rz[model] if r["d_new"] >= 30 and r["N"] >= 4 and not (0 < resonant_target_dim(r, "J") < 50)]
        D = np.array([current_edge_distance(r["s_new"]) for r in rs])
        jit = (np.array([r["N"] for r in rs]) - 5.5) * 0.004 + (0.006 if model == "generic" else -0.006)
        ax.scatter(D + jit, [r["kappa_max"] for r in rs], s=12, marker=mk, color=col, alpha=0.6,
                   label=f"max over sector ({MLAB[model]})")
        ax.scatter(D + jit, [r["kappa_mean"] for r in rs], s=12, marker=mk, facecolors="none", edgecolors=col, alpha=0.6,
                   label=f"mean over sector ({MLAB[model]})")
    ax.set_yscale("log")
    ax.set_xticks([0, 0.1, 0.2, 0.3])
    ax.set_xlabel(r"edge distance of the currents, $\min|s'\pm\frac{1}{5}|$")
    ax.set_ylabel(r"$\kappa$")
    ax.set_title(r"sectors with dim $B\geq30$, $N'=5$-$8$", fontsize=9)
    ax.legend(fontsize=6.5, loc="upper right")
    fig.tight_layout()
    return fig


def fig_edgegap(plt, EG):
    """Lowest levels of the chord-gapless family in the largest window-edge sector vs the family dimension."""
    fig, ax = plt.subplots(1, 1, figsize=(4.6, 3.8))
    for model, mk, col in (("sym", "o", "tab:blue"), ("generic", "^", "tab:green")):
        rs = EG[model]
        for N in sorted({r["N"] for r in rs}):
            big = max(r["dim"] for r in rs if r["N"] == N)
            sel = [r for r in rs if r["N"] == N and r["dim"] == big and r["low_rel"]]
            for j, r in enumerate(sel):
                D = r["family_dim"] * (1.03 if model == "generic" else 0.97)
                ax.scatter([D] * len(r["low_rel"][1:6]), r["low_rel"][1:6], s=6, marker=mk, color=col, alpha=0.25)
                ax.scatter(D, r["low_rel"][0], s=16, marker=mk, color=col, alpha=0.8,
                           label=MLAB[model] if (N == 7 and j == 0) else None)
    D = np.logspace(0.3, 4.2, 50)
    ax.plot(D, 0.6 / D, "k:", lw=0.9, label=r"$\propto D^{-1}$")
    ax.plot(D, 4.0 / D ** 2, "k--", lw=0.9, label=r"$\propto D^{-2}$ (hard edge)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("dimension $D$ of the gapless family (largest edge sector)")
    ax.set_ylabel("lowest levels / sector mean energy")
    ax.set_title("window edge: closes up to $N=7$, stalls at $N=9$", fontsize=9)
    ax.legend(fontsize=7)
    fig.tight_layout()
    return fig


def fig_tower(plt, T, Pr, TW=None):
    fig, axs = plt.subplots(1, 3 if TW else 2, figsize=(15 if TW else 10.5, 3.9))
    ax = axs[0]
    for lst, col, lab in ((T["other"], "tab:blue", r"target $|s|\neq\frac{1}{2}$"), (T["edge"], "tab:red", r"target $s=\pm\frac{1}{2}$ (gapless)")):
        ax.scatter([r["bh"] for r in lst], [r["tower"] for r in lst], s=16, color=col, alpha=0.75, label=lab)
    lo, hi = 3e-3, 3e2
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"black holes: $\langle h|O^\dagger G O|h\rangle$ averaged over $B$")
    ax.set_ylabel(r"tower: $\langle t|O^\dagger G O|t\rangle$ (one term)")
    ax.set_title(f"individual emissions; total $J$ cancels to {T['resid_max']:.0e}")
    ax.legend(fontsize=7)
    ax = axs[1]
    for model, mk, col in (("sym", "o", "tab:blue"), ("generic", "^", "tab:green")):
        rs = [r for r in Pr[model] if r.get("d_fort", 0) >= 10 and r["A"] == r["B"]]
        xs = [r["N"] + 1 + (r["A"] / (r["N"] + 1) - 0.5) * 0.6 for r in rs]
        ax.errorbar(xs, [r["c_state_q"][1] for r in rs], yerr=[[r["c_state_q"][1] - r["c_state_q"][0] for r in rs],
                    [r["c_state_q"][2] - r["c_state_q"][1] for r in rs]], fmt=mk, color=col, ms=3, lw=0.8,
                    label=f"fortuitous, {MLAB[model]} (median, 5-95%)")
    tw = [r for r in Pr["sym"] if "tower_c" in r and 0 < r["A"] < r["N"] + 1]
    ax.scatter([r["N"] + 1 + (r["A"] / (r["N"] + 1) - 0.5) * 0.6 for r in tw], [r["tower_c"] for r in tw], marker="*", s=40,
               color="tab:red", label=r"tower: $n(N'-n)/N'^2$")
    ax.axhline(0, color="0.5", lw=0.6)
    ax.set_xlabel(r"$N'$ (diagonal sectors $(n,n)$ spread by $n/N'$)")
    ax.set_ylabel(r"$c=\langle n_{\psi_e}n_{\chi_f}\rangle-\langle n_{\psi_e}\rangle\langle n_{\chi_f}\rangle$")
    ax.set_title("absorption of the new pair")
    ax.legend(fontsize=7)
    if TW:
        S, G = TW
        ax = axs[2]
        key = lambda r: (r["N"], r["seed"], r["n"])
        gd = {key(r): r for r in G}
        pts = [(gd[key(r)]["c2_mean"], r["c2_mean"], abs(abs(r["s_target"]) - 0.5) < 1e-9) for r in S if key(r) in gd]
        x, y, e = (np.array(v) for v in zip(*pts))
        ax.scatter(x[~e], y[~e], s=14, color="tab:blue", alpha=0.7, label="black holes, target off the edge")
        ax.scatter(x[e], y[e], s=14, color="tab:red", alpha=0.7, label=r"black holes, target $s=\pm\frac{1}{2}$")
        ax.scatter([1.0], [1.0], marker="*", s=90, color="k", label=r"tower: $\cos^2=1$ (all $N,n$)")
        ax.plot([1e-3, 1.2], [1e-3, 1.2], "k--", lw=0.8)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(8e-3, 1.5)
        ax.set_ylim(8e-3, 1.5)
        ax.set_xlabel(r"generic model: $\langle\cos^2\rangle(EB_{00},KB_{11})$")
        ax.set_ylabel(r"symmetrized model (tower removed)")
        ax.set_title("twisted chord: no coherence for black holes")
        ax.legend(fontsize=6.5, loc="upper left")
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-figs", action="store_true")
    a = ap.parse_args()
    summary, W, B, Rz, Pr = {}, {}, {}, {}, {}
    for model in MODELS:
        summary[f"window_{model}"], W[model] = check_window(model)
        summary[f"blocks_{model}"], B[model] = check_blocks(model)
        summary[f"resonance_{model}"], Rz[model] = check_resonance(model)
        summary[f"pairing_{model}"], Pr[model] = check_pairing(model)
    summary["tower"], _ = check_tower()
    summary["twisted"], TWs, TWg = check_twisted()
    EG = {m: load("edgegap", m) for m in MODELS}
    summary["edgegap"] = {m: sorted({(r["N"], r["seed"], r["a"], r["b"], r["dim"], r["family_dim"], round(r["gap_rel"], 10))
                                     for r in EG[m] if r["gap_rel"]}) for m in MODELS}
    T = summary["tower"]
    for model in MODELS:
        print(f"\n===== {model}")
        for Np, v in summary[f"window_{model}"]["by_Np"].items():
            print(f"[window] N'={Np}: BPS sectors {v['n_bps_sectors']}/{v['n_sectors']}; outside window (not edge): {v['outside_window']}; "
                  f"frac/D mean,min,max={v['frac_over_D']} central={v['frac_over_D_central']}; edge gap/m1 range={v['edge_gap_rel']}, "
                  f"bulk median/min={v['bulk_gap_rel_median']}/{v['bulk_gap_rel_min']}" + (f"; boundary={v.get('boundary')}" if Np >= 9 else ""))
        for N, v in summary[f"blocks_{model}"].items():
            print(f"[blocks] N={N} lam={v['lam']:.2f}: point mass gm/rms/n={v['point_mass']}, continuum far={v['cont_far']}, "
                  f"near={v['cont_near']}, large-edge-target T/c0 median {v['edge_T_over_c0_median']} vs off-edge {v['other_T_over_c0_median']}")
        for Np, v in summary[f"resonance_{model}"].items():
            print(f"[resonance] N'={Np} seeds {v['seeds']}: J-resonant (large target) {v['J_resonant_large_target']}; "
                  f"J-resonant (small target) {v['J_resonant_small_target']}; z-resonant {v['z_resonant']}; "
                  f"non-resonant {v['nonresonant']}; tower bare fraction min {v['tower_bare_fraction_min']}")
        for Np, v in summary[f"pairing_{model}"].items():
            print(f"[pairing] N'={Np}: fortuitous c mean range {v['c_mean_range']} median {v['c_mean_median']}, max |<n_psi>-A/N'| "
                  f"{v['n_psi_dev_max']}; tower c {v['tower_c']}")
    for k, v in summary["twisted"].items():
        print(f"[twisted] {k} s_t={v['s_target']:+.1f}: <cos^2> sym {v['sym'][0]:.4f}+-{v['sym'][1]:.4f}, generic "
              f"{v['generic'][0] if v['generic'] else float('nan'):.4f}; tower cos^2 min {v['tower_cos2_min']:.12f}")
    print(f"\n[tower] {T['n']} records, max residual {T['resid_max']:.2e}")
    for r in T["edge"]:
        print(f"   edge target: N={r['N']} seed={r['seed']} n={r['n']} {r['current']} s={r['s_target']:+.1f} gap/m1={r['gap_rel']:.2e}: "
              f"tower {r['tower']:.3f} vs black-hole mean {r['bh']:.3f}")
    with open(os.path.join(DATA, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1, default=lambda o: o.item() if isinstance(o, np.generic) else str(o))
    if a.no_figs:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.labelsize": 10, "legend.fontsize": 8, "font.family": "DejaVu Sans"})
    os.makedirs(FIGS, exist_ok=True)
    for name, fig in (("fig_window", fig_window(plt, W)), ("fig_dictionary", fig_dictionary(plt, B, W)),
                      ("fig_resonance", fig_resonance(plt, Rz)), ("fig_tower", fig_tower(plt, T, Pr, (TWs, TWg) if TWs else None)), ("fig_edgegap", fig_edgegap(plt, EG))):
        for ext in ("pdf", "png"):
            fig.savefig(os.path.join(FIGS, f"{name}.{ext}"), dpi=200 if ext == "png" else None)
        plt.close(fig)
        print("wrote", name)


if __name__ == "__main__":
    main()
