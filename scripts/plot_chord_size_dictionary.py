#!/usr/bin/env python3
r"""
Figure for research/notes/chord_size_dictionary_2026-10-02.md (chord number vs operator size of the BPS projector).

Panels (all numbers recomputed from raw data):
 (a) chord-number generating function E tau(P(omega) P(omega_x)) vs x^2 (SYK, angle-rotated OU pairs) against BLY's
     sum_n P_n x^{2n} at lambda = 2 p^2 / N' (BLY eqs. 3.10, 4.5-4.6 with q^Delta -> x);
 (b) exact chord-1 weight ratio pi_1/pi_0 relative to BLY's 2q/(1-q), as p -> infinity at fixed lambda;
 (c) size composition J_nk / pi_n of each chord sector at N' = 12 (free fits of exact size overlaps);
 (d) diagonal share J_{n,2n} / pi_n vs N' (n = 1 exact; N' = 14 from top-size projections and Krylov pi_n).

Inputs: results/data/chord_chaos_2026-10-01/{exact, exact_K16, topsize}/ ; outputs results/figures/chord_size_dictionary_2026-10-02/.
Usage: .venv/bin/python scripts/plot_chord_size_dictionary.py
"""
import glob
import json
import os
import sys
from math import comb, exp

import matplotlib
import matplotlib.ticker

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from analyze_chord_chaos import krylov_kappa, load, lsq, pool_lags  # noqa: E402
from chord_chaos import pi1_ratio_exact  # noqa: E402
from dssyk_n2 import hh_length_distribution  # noqa: E402

DATA = os.path.join(ROOT, "results/data/chord_chaos_2026-10-01")
OUT = os.path.join(ROOT, "results/figures/chord_size_dictionary_2026-10-02")

# validated palette (dataviz reference instance): categorical slots 1-3, ordinal blue ramp for size k
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"
SIZE_RAMP = {0: "#86b6ef", 2: "#3987e5", 4: "#1c5cab", 6: "#0d366b"}
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"


def style(ax):
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#c3c2b7")
    ax.tick_params(colors=INK2, labelsize=9)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)


def free_joint(N):
    """J_nk (rows k, cols n) by free polynomial fits of exact size overlaps; returns J, t-grid, pooled G."""
    rows, K, G, Gt, Ct, Cf, counts, Sk = load(os.path.join(DATA, "exact", f"chord_chaos_p3_N{N}.jsonl"))
    t = np.cos(np.pi * np.arange(K // 2 + 1) / K) ** 2
    rmax = Sk.shape[-1] - 1
    J = np.zeros((rmax + 1, 7))
    for k in range(0, rmax + 1, 2):
        J[k] = lsq(t, pool_lags(Sk[:, :, k], K).mean(0), list(range(7)))
    Gp = pool_lags(G, K)
    return J, t, Gp.mean(0), Gp.std(0, ddof=1) / np.sqrt(len(Gp)), rows[0]["a"]


def main():
    os.makedirs(OUT, exist_ok=True)
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 10, "axes.titlesize": 10.5})
    fig, axs = plt.subplots(2, 2, figsize=(10.2, 7.8))
    summary = {}

    # (a) OU generating function vs BLY
    ax = axs[0, 0]
    tt = np.linspace(0, 1, 200)
    for N, col, mk in ((8, C1, "o"), (10, C2, "s"), (12, C3, "^")):
        J, t, g, ge, a = free_joint(N)
        Pn = hh_length_distribution(18 / N, 0.0)
        bly = np.array([sum(Pn[n] * x ** n for n in range(len(Pn))) for x in tt])
        ax.plot(tt, bly, color=col, lw=1.5, ls="--")
        ax.plot(tt, a + (1 - a) * (bly - Pn[0]) / (1 - Pn[0]), color=col, lw=1.0, ls=":")
        ax.errorbar(t, g, yerr=ge, fmt=mk, ms=6, color=col, mec="white", mew=0.8, elinewidth=1, zorder=3)
        ax.text(-0.025, g[-1], f"$N'={N}$", ha="right", va="center", color=INK2, fontsize=9)
        summary[f"G_N{N}"] = dict(t=t.tolist(), G=g.tolist(), G_err=ge.tolist(), a=a)
    ax.set_xlabel(r"$x^2$  (coupling correlation squared)")
    ax.set_ylabel(r"$E\,\tau\,(P(\omega)P(\omega_x)) = \sum_n \pi_n x^{2n}$")
    ax.set_xlim(-0.2, 1.03)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_title("(a) Chord-number generating function: SYK (points),\n"
                 "BLY $\\sum_n P_n x^{2n}$ (dashed), BLY with $P_0\\to a$ (dotted)", loc="left", color=INK)
    style(ax)

    # (b) double-scaled limit of the exact chord-1 weight
    ax = axs[0, 1]
    ps = np.array([3, 5, 7, 9, 13, 17, 25, 35])
    for lam, col, mk in ((1.0, C1, "o"), (2.0, C2, "s")):
        dev, pp = [], []
        for p in ps:
            N = round(2 * p * p / lam)
            N += N % 2
            lam_e = 2 * p * p / N
            q = exp(-lam_e)
            dev.append(float(pi1_ratio_exact(N, int(p))) / (2 * q / (1 - q)) - 1)
            pp.append(p)
        dev = np.array(dev)
        ax.loglog(pp, dev, mk + "-", color=col, lw=1.5, ms=6, mec="white", mew=0.8, label=fr"$\lambda={lam:g}$")
        ax.loglog(ps, dev[-1] * (ps[-1] / ps) ** 2, ":", color=MUTED, lw=1)
        summary[f"pi1_dev_lam{lam:g}"] = dict(p=pp, rel_dev=dev.tolist())
    ax.annotate(r"$\propto p^{-2}$", (9, 0.012), color=MUTED, fontsize=9)
    ax.set_xticks(ps, [str(x) for x in ps])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.annotate("our model ($p=3$)", (3, 0.0629), xytext=(4.2, 0.12), color=INK2, fontsize=9,
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
    ax.set_xlabel(r"$p$  (fermions per supercharge term), $N'=2p^2/\lambda$")
    ax.set_ylabel(r"$(\pi_1/\pi_0)_{\rm exact}\,/\,[2q/(1-q)] - 1$")
    ax.set_title("(b) Exact chord-1 weight $\\to$ BLY in double scaling", loc="left", color=INK)
    ax.legend(frameon=False, fontsize=9, labelcolor=INK2)
    style(ax)

    # (c) size composition of each chord sector at N' = 12
    ax = axs[1, 0]
    J12, *_ = free_joint(12)
    ns = [1, 2, 3, 4]
    bottoms = np.zeros(len(ns))
    width = 0.62
    comp = {}
    for k in (2, 4, 6):
        frac = np.array([max(J12[k, n], 0) / max(J12[:, n].clip(min=0).sum(), 1e-12) for n in ns])
        ax.bar(ns, frac, width, bottom=bottoms, color=SIZE_RAMP[k], edgecolor="white", linewidth=1.5,
               label=fr"size $k={k}$")
        for x, b, f in zip(ns, bottoms, frac):
            if f > 0.06:
                ax.text(x, b + f / 2, f"{f:.2f}", ha="center", va="center", fontsize=8.5,
                        color="white" if k >= 4 else INK)
        comp[k] = frac.tolist()
        bottoms += frac
    ax.set_xticks(ns, [f"$n={n}$\n$k\\leq{2*n}$" for n in ns])
    ax.set_ylim(0, 1.2)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_ylabel(r"$J_{nk}/\pi_n$  (size composition of chord sector $n$)")
    ax.set_title("(c) Chord sectors at $N'=12$: support $k\\leq 2n$ (exact),\n"
                 "top size dominant until the cap $k\\leq F=6$", loc="left", color=INK)
    ax.legend(frameon=False, fontsize=9, labelcolor=INK2, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.02))
    style(ax)
    summary["composition_N12"] = dict(n=ns, **{f"k{k}": v for k, v in comp.items()})

    # (d) diagonal share vs N'
    ax = axs[1, 1]
    Ns = [8, 10, 12]
    s2, s3 = [], []
    for N in Ns:
        J, *_ = free_joint(N)
        s2.append(J[4, 2] / J[:, 2].sum())
        s3.append(J[6, 3] / J[:, 3].sum() if J.shape[0] > 6 else np.nan)
    # N' = 14: J_{2,4}, J_{3,6} from top-size projections (lower bounds), pi_n from Krylov (pi_2, pi_3 ~ 1.03-1.04
    # kappa_2 and ~0.93-0.95 kappa_3 at N' = 10, 12) -> central value and a systematic band
    rows14 = [json.loads(l) for l in open(os.path.join(DATA, "topsize", "topsize_p3_N14.jsonl"))]
    r14 = np.array([r["r"] for r in rows14]); d14 = rows14[0]["d"]
    m = r14.mean(0)
    J24_14, J36_14 = m[0] ** 2 / (d14 * m[1]), m[2] ** 2 / (d14 * m[3])
    specs = sorted(glob.glob(os.path.join(DATA, "exact*", "spectra_p3_N14_s*.npy")))
    kap14, _, _ = krylov_kappa(specs, 2 * 3 ** 6 / comb(14, 7), 5)
    pi2_lo, pi2_hi = kap14[2], 1.05 * kap14[2]
    pi3_lo, pi3_hi = 0.93 * kap14[3], kap14[3]
    s2_14 = (J24_14 / pi2_hi, J24_14 / pi2_lo)
    s3_14 = (J36_14 / pi3_hi, J36_14 / pi3_lo)
    ax.plot([8, 14], [1, 1], color=C3, lw=1.5)
    ax.annotate(r"$n=1$: exact ($J_{1k}=\pi_1\delta_{k2}$)", (8.1, 1.015), color=INK2, fontsize=9)
    ax.plot(Ns, s2, "o-", color=C1, lw=1.5, ms=7, mec="white", mew=0.8)
    ax.errorbar([14], [np.mean(s2_14)], yerr=[[np.mean(s2_14) - s2_14[0]], [s2_14[1] - np.mean(s2_14)]],
                fmt="o", mfc="white", mec=C1, color=C1, ms=7, elinewidth=1.2, capsize=3)
    ax.plot([Ns[-1], 14], [s2[-1], np.mean(s2_14)], ":", color=C1, lw=1.2)
    ax.annotate(r"$n=2$ ($k=4$)", (10.15, s2[1] - 0.07), color=INK2, fontsize=9)
    ax.plot([12], [s3[-1]], "s", color=C2, ms=7, mec="white", mew=0.8)
    ax.errorbar([14], [np.mean(s3_14)], yerr=[[np.mean(s3_14) - s3_14[0]], [s3_14[1] - np.mean(s3_14)]],
                fmt="s", mfc="white", mec=C2, color=C2, ms=7, elinewidth=1.2, capsize=3)
    ax.plot([12, 14], [s3[-1], np.mean(s3_14)], ":", color=C2, lw=1.2)
    ax.annotate(r"$n=3$ ($k=6$)" + "\n" + r"$N'=12$: $2n=F$ (size cap)", (10.6, 0.13), color=INK2, fontsize=9)
    ax.set_xticks([8, 10, 12, 14])
    ax.set_xlabel(r"$N'$  (complex fermions; half filling $F=N'/2$)")
    ax.set_ylabel(r"$J_{n,2n}/\pi_n$  (share of chord sector $n$ at size $k=2n$)")
    ax.set_ylim(0, 1.08)
    ax.set_title("(d) Diagonal $k=2n$ becomes dominant with $N'$\n(open: $N'=14$ from top-size projections)",
                 loc="left", color=INK)
    style(ax)
    summary["diag_share"] = dict(N=Ns + [14], n2=s2 + [list(s2_14)], n3=s3[:2] + [s3[2], list(s3_14)],
                                 J24_N14_proj=J24_14, J36_N14_proj=J36_14, kappa_N14=kap14.tolist())

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"fig_chord_size_dictionary.{ext}"), dpi=200)
    with open(os.path.join(OUT, "fig_chord_size_dictionary_data.json"), "w") as fh:
        json.dump(summary, fh, indent=1, default=float)
    print(json.dumps({k: v for k, v in summary.items() if k.startswith("diag") or k.startswith("comp")}, default=float))


if __name__ == "__main__":
    main()
