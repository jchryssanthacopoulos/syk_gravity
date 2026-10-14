#!/usr/bin/env python3
r"""
Figures for the chord computation of the first-order dressing spectrum (progress report 2026-10-14, docs/derivations.md D7).

  fig_weight         BPS -> continuum spectral weight of the pair operator: closed form (Al-Salam-Chihara) vs finite N
  fig_t1             mean first-order dressing cost t1 vs N': exact (mean +- SE) vs 2phi1 closed form
  fig_W1             W^(1) vs N': exact first order, chord prediction, planar law with exact inputs, decoder; error budget
  fig_sigma2         distribution of first-order dressing costs: exact vs planar (Kronecker) vs chord
  fig_heterogeneity  BPS heterogeneity of the averaged dressing cost (rsd vs N', tau dependence, free estimate)

Inputs: results/data/chord_dressing_2026-10-14/{decomp,heterogeneity}_N*.jsonl, chord_predictions.json,
        results/data/uplift_2026-10-12/walls_N*.jsonl (exact first order, 13 realizations) and the decoder W from the
        2026-10-13 report table (hard-coded below with its source).
Usage: .venv/bin/python scripts/plot_chord_dressing.py [--only fig_W1 ...]
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
DATA = os.path.join(ROOT, "results/data/chord_dressing_2026-10-14")
WALLS = os.path.join(ROOT, "results/data/uplift_2026-10-12")
FIGS = os.path.join(ROOT, "results/figures/chord_dressing_2026-10-14")

# decoder heterogeneity W^dec = T2^dec - T4^dec, all-site OU campaign (progress report 2026-10-13, Table 5);
# N' = 16: exact decoder on the two pre-registered realizations (W^(1) + residual +0.0030, +0.0042; same report, Table 6)
W_DEC = {8: 0.0396, 9: 0.0653, 10: 0.0666, 11: 0.0891, 12: 0.0897, 13: 0.1083, 14: 0.1085, 15: 0.1230}


def load_jsonl(path):
    return [json.loads(l) for l in open(path)] if os.path.exists(path) else []


def walls_W1(N):
    rows = {}
    for r in load_jsonl(os.path.join(WALLS, f"walls_N{N}.jsonl")):
        rows[r["seed"]] = {**rows.get(r["seed"], {}), **r}
    W, T = [], []
    for r in rows.values():
        y = (2 * np.array(r["first_lam"]) - 1) ** 2
        W.append(np.mean(y * (1 - y)))
        T.append(np.mean(y))
    return np.array(W), np.array(T)


def decomp(N, tag=""):
    return load_jsonl(os.path.join(DATA, f"decomp_N{N}{tag}.jsonl"))


def ms(v):
    v = np.asarray(v, float)
    return v.mean(), (v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else np.nan)


def fig_weight(plt):
    from chord_dressing import DressingFamily, chord_sectors, eth_function, old_model_lean
    from dssyk_n2 import bps_continuum_weight, chord_energy
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    for ax, Np in zip(axes, (14, 13)):
        P = chord_sectors(Np, "h")
        Es, ws = [], []
        for seed in (9000, 9001, 9002, 9003):
            Fm = DressingFamily(old_model_lean(Np, seed), "h")
            E, Vnn = eth_function(Fm)
            Es.append(E)
            ws.append(Vnn / (Fm.b * P["c0"]))          # weight of each inner level in the normalized two-point function
        E, w = np.concatenate(Es), np.concatenate(ws) / 4
        edges = np.linspace(0, np.percentile(E, 99.8) * 1.02, 45)
        h, _ = np.histogram(E, bins=edges, weights=w)
        ax.stairs(h / np.diff(edges), edges, color="k", lw=1.3, label=f"exact, $N'={Np}$ (4 realizations)")
        phi = np.linspace(1e-4, np.pi - 1e-4, 4000)
        Ech = P["kappa"] * chord_energy(P["lam"], P["s_in"], phi)
        f = bps_continuum_weight(P["lam"], P["s_out"], P["delta"], phi)
        dEdphi = P["kappa"] * 2 * P["q"] ** -0.5 * np.sin(phi)
        ax.plot(Ech, f / dEdphi, color="tab:red", lw=1.6, label="chord (Al-Salam–Chihara), no fit")
        ax.axvline(P["kappa"] * chord_energy(P["lam"], P["s_in"], 0.0), color="tab:red", ls=":", lw=1)
        ax.set_xlabel(r"inner energy $E$ (charge $s_{\rm in}=%.3g$)" % P["s_in"])
        ax.set_ylabel(r"spectral weight $f(E)$")
        ax.set_title(rf"$N'={Np}$: $\lambda={P['lam']:.3f}$, $s_{{\rm out}}={P['s_out']:.3g}$", fontsize=10)
        ax.legend(frameon=False, fontsize=8)
        ax.set_xlim(0, edges[-1])
    fig.tight_layout()
    return fig


def fig_t1(plt):
    from chord_dressing import chord_t1
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    Ns = list(range(10, 18))
    for fam, col, mk in (("h", "tab:blue", "o"), ("g", "tab:orange", "s")):
        xs, m, s = [], [], []
        for N in Ns:
            R = decomp(N)
            if not R or fam not in R[0]:
                continue
            mm, ss = ms([r[fam]["t1"] for r in R])
            xs.append(N)
            m.append(mm)
            s.append(ss)
        ax.errorbar(xs, m, yerr=s, fmt=mk, color=col, ms=5, capsize=2, label=f"exact, {fam} family")
        Nc = np.arange(10, 18)
        ax.plot(Nc, [chord_t1(int(n), fam) for n in Nc], ls="--", color=col, lw=1.2, marker=mk, mfc="none", ms=5,
                label=f"chord $_2\\phi_1$, {fam} family")
    ax.set_xlabel(r"$N'$")
    ax.set_ylabel(r"$t_1=\langle\sigma^2\rangle=c_0\int_0^\infty G(\tau)\,d\tau$")
    ax.legend(frameon=False, fontsize=8, ncol=2)
    fig.tight_layout()
    return fig


def fig_W1(plt):
    pred = json.load(open(os.path.join(DATA, "chord_predictions.json")))["predictions"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.0), gridspec_kw=dict(width_ratios=[1.35, 1]))
    ax = axes[0]
    Ns = range(10, 18)
    ex = {N: walls_W1(N)[0] for N in Ns}
    xe = [N for N in Ns if len(ex[N])]
    ax.errorbar(xe, [ex[N].mean() for N in xe], yerr=[ex[N].std(ddof=1) / np.sqrt(len(ex[N])) if len(ex[N]) > 1 else 0 for N in xe],
                fmt="o", color="k", ms=5, capsize=2, label=r"exact first order $W^{(1)}$")
    x17 = [r for r in decomp(17)]
    if x17:
        m, s = ms([r["own"]["W"] for r in x17])
        ax.errorbar([17], [m], yerr=[s], fmt="o", color="k", mfc="white", ms=6, capsize=2, label=r"$N'=17$ (pre-registered test)")
    pre = os.path.join(DATA, "prereg_N17.json")
    if os.path.exists(pre):
        P1 = [p for p in json.load(open(pre))["predictions"] if p["id"] == "P1"][0]
        lo, hi = P1["pass_interval"]
        ax.errorbar([17], [P1["central"]], yerr=[[P1["central"] - lo], [hi - P1["central"]]], fmt="_", color="0.3", ms=12,
                    capsize=4, lw=1.4, label=r"$N'=17$: pre-registered $W^{(1)}$ interval")
    xs = sorted(int(k) for k in pred)
    xs_s = [x for x in xs if x <= 20]
    ax.plot(xs_s, [pred[str(x)]["W"] for x in xs_s], "--", color="tab:red", marker="^", ms=5, lw=1.2,
            label="chord prediction (no fit)")
    kr = [(N, ms([r["kron"]["W"] for r in decomp(N)])[0]) for N in range(10, 18) if decomp(N) and "kron" in decomp(N)[0]]
    ax.plot([a for a, _ in kr], [b for _, b in kr], ":", color="tab:blue", marker="s", ms=4, lw=1,
            label="planar law, exact inputs")
    xd = sorted(W_DEC)
    xd = [x for x in xd if x >= 10]
    ax.plot(xd, [W_DEC[x] for x in xd], "D", color="tab:green", ms=5, label=r"decoder $W^{\rm dec}$")
    ax.axvspan(16.6, 17.4, color="0.92", zorder=0)
    ax.set_xlabel(r"$N'$")
    ax.set_ylabel(r"$W=\langle y(1-y)\rangle$")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax = axes[1]
    rows = []
    for N in range(10, 18):
        R = decomp(N)
        if not R or "kron" not in R[0]:
            continue
        S = np.mean([r["sampled"]["W"] for r in R])
        rows.append((N, (np.mean([r["kron"]["W"] for r in R]) / S - 1) * 100,
                     (np.mean([r["chord"]["W"] for r in R]) / np.mean([r["wishart"]["W"] for r in R]) - 1) * 100,
                     (np.mean([r["wishart"]["W"] for r in R]) / np.mean([r["kron"]["W"] for r in R]) - 1) * 100,
                     (np.mean([r["chord"]["W"] for r in R]) / S - 1) * 100))
    if rows:
        x = np.array([r[0] for r in rows])
        ax.bar(x - 0.27, [r[1] for r in rows], width=0.27, color="tab:blue", label="planar law (exact inputs)")
        ax.bar(x, [r[2] for r in rows], width=0.27, color="tab:red", label="chord inputs")
        ax.bar(x + 0.27, [r[3] for r in rows], width=0.27, color="0.6", label="BPS factorization")
        ax.plot(x, [r[4] for r in rows], "k^-", ms=5, lw=1, label="total (chord vs exact)")
        ax.axhline(0, color="k", lw=0.5)
        ax.set_xlabel(r"$N'$")
        ax.set_ylabel(r"relative error in $W^{(1)}$ (%)")
        ax.set_ylim(-2.5, 12)
        ax.legend(frameon=False, fontsize=7.5, loc="upper left", ncol=2)
    fig.tight_layout()
    return fig


def fig_sigma2(plt):
    from chord_dressing import DressingFamily, analyze_family, old_model_lean, wy
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    for ax, Np in zip(axes, (14, 15)):
        om = old_model_lean(Np, 9000)
        rng = np.random.default_rng(1)
        Fm = DressingFamily(om, "h")
        _, sp = analyze_family(Fm, om["eta"], 8 if Np <= 14 else 4, rng, draws=8 if Np <= 14 else 4,
                               third_moment=False, Np=Np)
        bins = np.logspace(-3.5, 1.7, 45)
        for key, col, lab in (("sampled", "k", "exact first order"), ("kron", "tab:blue", "planar law, exact inputs"),
                              ("chord", "tab:red", "chord prediction")):
            x = np.clip(sp[key], 1e-12, None)
            h, _ = np.histogram(x, bins=bins)
            ax.stairs(h / len(x) / np.diff(np.log10(bins)), bins, color=col, lw=1.3, label=lab)
        s = np.logspace(-3.5, 1.7, 300)
        ax2 = ax.twinx()
        ax2.plot(s, wy(s)[0], color="0.6", lw=0.8, ls="--")
        ax2.set_ylim(0, 0.6)
        ax2.set_yticks([])
        ax.set_xscale("log")
        ax.set_xlabel(r"first-order dressing cost $\sigma^2$ ($h$ family)")
        ax.set_ylabel(r"density per decade")
        ax.set_title(rf"$N'={Np}$ (dashed grey: weight $w(\sigma^2)=y(1-y)$)", fontsize=10)
        ax.legend(frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout()
    return fig


def fig_heterogeneity(plt):
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    Ns = list(range(10, 18))
    for fam, col, mk in (("h", "tab:blue", "o"), ("g", "tab:orange", "s")):
        x, m, f, cb = [], [], [], []
        for N in Ns:
            R = load_jsonl(os.path.join(DATA, f"heterogeneity_N{N}.jsonl"))
            if not R or fam not in R[0]:
                continue
            x.append(N)
            m.append(np.mean([r[fam]["rsd"] for r in R]))
            f.append(np.mean([r[fam]["free"]["rsd_free_all"] for r in R]))
            ir_o = np.mean([r[fam]["bins"]["rsd_bin"][0] for r in R])
            ir_f = np.mean([r[fam]["free"]["rsd_free_bins"][0] for r in R])
            cb.append(max(0.0, 1 - (ir_f / ir_o) ** 2))
        if fam == "g":
            keep = [i for i, n in enumerate(x) if n % 2 == 1]           # g = h by Hodge duality at even N'
            x, m, f, cb = [x[i] for i in keep], [m[i] for i in keep], [f[i] for i in keep], [cb[i] for i in keep]
        axes[0].plot(x, m, marker=mk, color=col, lw=1.2, label=rf"observed, {fam}")
        axes[0].plot(x, f, marker=mk, color=col, lw=1, ls=":", mfc="none", label=rf"free (incoherent) estimate, {fam}")
        axes[2].plot(x, cb, marker=mk, color=col, lw=1.2, label=rf"{fam} family")
    axes[0].set_xlabel(r"$N'$")
    axes[0].set_ylabel(r"rel. spread of $\bar{\mathcal{D}}=E_\eta\mathcal{D}$ over BPS states")
    axes[0].set_ylim(0, 0.6)
    axes[0].legend(frameon=False, fontsize=7.2, loc="upper left", ncol=2)
    for N, ls in ((11, "-"), (13, "--"), (15, ":"), (12, "-"), (14, "--"), (16, ":")):
        R = load_jsonl(os.path.join(DATA, f"heterogeneity_N{N}.jsonl"))
        if not R:
            continue
        taus = [float(t) for t in R[0]["h"]["tau"]]
        v = np.mean([[r["h"]["tau"][str(t)] for t in R[0]["h"]["tau"]] for r in R], axis=0)
        axes[1].plot(taus, v, ls=ls, color="tab:blue" if N % 2 else "tab:purple", marker="o", ms=3, label=rf"$N'={N}$")
    axes[1].set_xscale("log")
    axes[1].set_xlabel(r"$\tau E_{\rm med}$")
    axes[1].set_ylabel(r"rel. spread of $K(\tau)=P\,E_\eta[E^\dagger e^{-\tau H_0}E]\,P$")
    axes[1].legend(frameon=False, fontsize=7.5, ncol=2)
    axes[2].set_xlabel(r"$N'$")
    axes[2].set_ylabel(r"coherent share, lowest 5% of levels")
    axes[2].set_ylim(0, 1)
    axes[2].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    return fig


def tables():
    """LaTeX rows for the report tables (t1, error budget, heterogeneity)."""
    from chord_dressing import chord_sectors
    from dssyk_n2 import chord_energy
    pred = json.load(open(os.path.join(DATA, "chord_predictions.json")))["predictions"]
    print("% TABLE_T1")
    for N in range(10, 18):
        R = decomp(N)
        if not R:
            continue
        for fam in ("h", "g"):
            if fam == "g" and N % 2 == 0:
                continue                                           # g = h by Hodge duality at even N'
            m, s = ms([r[fam]["t1"] for r in R])
            P = chord_sectors(N, fam)
            gap = P["kappa"] * chord_energy(P["lam"], P["s_in"], 0.0)
            lab = f"{N}, {fam}" + (" (= $g$)" if fam == "h" and N % 2 == 0 else "")
            mtxt = f"${m:.3f}\\pm{s:.3f}$" if np.isfinite(s) else f"${m:.3f}$"
            print(f"{lab} & {len(R)} & {mtxt} & {R[0][fam]['t1_chord']:.3f} & {R[0][fam]['t1_chord'] / m:.2f} & "
                  f"{np.mean([r[fam]['diag_frac_S'] for r in R]):.4f} & {np.mean([r[fam]['Emin'] for r in R]):.1f} / {gap:.1f}\\\\")
    print("% TABLE_BUDGET")
    for N in range(10, 18):
        R = decomp(N)
        W, T = walls_W1(N)
        if not R:
            continue
        f = lambda k, q="W": np.mean([r[k][q] for r in R])           # noqa: E731
        wex = W.mean() if len(W) else f("own")
        tex = T.mean() if len(T) else f("own", "T2")
        wd = W_DEC.get(N, None)
        print(f"{N} & {wex:.4f} & {f('kron'):.4f} & {f('wishart'):.4f} & {f('chord'):.4f} & {f('chord') / wex:.3f} & "
              f"{tex:.3f} / {f('chord', 'T2'):.3f} & {wd if wd else '--'}\\\\")
    print("% HETERO")
    for N in range(10, 18):
        R = load_jsonl(os.path.join(DATA, f"heterogeneity_N{N}.jsonl"))
        if not R:
            continue
        for fam in ("h", "g"):
            if fam == "g" and N % 2 == 0:
                continue
            obs = np.mean([r[fam]["rsd"] for r in R])
            fr = np.mean([r[fam]["free"]["rsd_free_all"] for r in R])
            ir_o = np.mean([r[fam]["bins"]["rsd_bin"][0] for r in R])
            ir_f = np.mean([r[fam]["free"]["rsd_free_bins"][0] for r in R])
            vs = np.mean([r[fam]["bins"]["var_share"][0] for r in R])
            cb = np.mean([r[fam]["probes"]["Bc"] for r in R])
            rb = np.mean([r[fam]["probes"]["rsd_Bc"] for r in R])
            ne = np.mean([r[fam]["free"]["n_eff_levels"] for r in R])
            print(f"{N}, {fam} & {obs:.3f} & {fr:.3f} & {ir_o:.2f} / {ir_f:.2f} & {max(0.0, 1 - (ir_f / ir_o) ** 2):.2f} & {vs:.2f} & "
                  f"{ne:.0f} & {cb:+.2f} & {rb:.3f}\\\\")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--out", default=FIGS)
    ap.add_argument("--tables", action="store_true", help="print LaTeX table rows and exit")
    a = ap.parse_args()
    if a.tables:
        tables()
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.labelsize": 10, "legend.fontsize": 8})
    os.makedirs(a.out, exist_ok=True)
    for name in ("fig_weight", "fig_t1", "fig_W1", "fig_sigma2", "fig_heterogeneity"):
        if a.only and name not in a.only:
            continue
        fig = globals()[name](plt)
        path = os.path.join(a.out, name + ".pdf")
        fig.savefig(path)
        print("wrote", path, flush=True)


if __name__ == "__main__":
    main()
