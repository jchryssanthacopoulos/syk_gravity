#!/usr/bin/env python3
r"""
Walls collapse, localization does not (progress report 2026-10-13; docs/derivations.md D6 §§4, 6).

(a) Fraction of decoder eigenvalues within 1e-3 of a wall (min(lambda, 1 - lambda) < 1e-3), exact vs first-order
    uplift from the old model, and the fraction of exact atoms; right axis: the dimension ratio
    r = e_{F+2} / b_F (exact forms of degree F+2 over old BPS states of degree F) that controls the hard edge.
(b) Localization L = T4_dec - T4_OU(x*) measured (OU campaign, D5 + N'=15) and predicted from the first-order
    heterogeneity, L_pred = T2_dec - H_first - T4_OU(x*).

Inputs: results/data/uplift_2026-10-12/walls_N*.jsonl, heterogeneity_analysis.json.
Usage: .venv/bin/python scripts/plot_uplift_walls.py
"""
import argparse
import json
import os
import sys

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from uplift import atom_bound  # noqa: E402

WALLS = os.path.join(ROOT, "results/data/uplift_2026-10-12")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Np", type=int, nargs="+", default=list(range(8, 17)))
    ap.add_argument("--out", default=os.path.join(ROOT, "results/figures/uplift_2026-10-12/fig_walls_vs_localization.pdf"))
    ap.add_argument("--sym", default="H", help="symbol for <y(1-y)> in legends (the 2026-10-13 report uses W)")
    a = ap.parse_args()
    rows = []
    for N in a.Np:
        R = [json.loads(l) for l in open(os.path.join(WALLS, f"walls_N{N}.jsonl"))]
        f1 = np.array([r["first"][3] for r in R])
        ex = np.array([r["exact"][3] for r in R if "exact" in r])
        at = np.array([r["exact_atoms"] / len(r["lam"]) for r in R if "exact" in r])
        rows.append(dict(N=N, first=f1, exact=ex, atoms=at, r=atom_bound(N)["ratio_h"]))
        print(N, f"first {f1.mean():.4f} [{f1.min():.4f},{f1.max():.4f}]  exact {ex.mean():.4f} [{ex.min():.4f},{ex.max():.4f}]"
                 f"  atoms {at.mean():.4f}  r {rows[-1]['r']:.3f}")
    het = {e["N"]: e for e in json.load(open(os.path.join(WALLS, "heterogeneity_analysis.json")))}

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.labelsize": 10, "legend.fontsize": 8.2})
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.7))
    N = np.array([r["N"] for r in rows])

    def band(key):
        m = np.array([r[key].mean() for r in rows])
        lo = np.clip(m - np.array([r[key].min() for r in rows]), 0, None)
        hi = np.clip(np.array([r[key].max() for r in rows]) - m, 0, None)
        return m, [lo, hi]
    m, err = band("exact")
    ax[0].errorbar(N, m, yerr=err, fmt="ko-", ms=5, capsize=2, label=r"exact decoder")
    m1, err1 = band("first")
    ax[0].errorbar(N, m1, yerr=err1, fmt="D--", color="tab:red", mfc="none", ms=5, capsize=2,
                   label="first-order uplift (old model)")
    ma, _ = band("atoms")
    pos = ma > 0
    ax[0].plot(N[pos], ma[pos], "s:", color="tab:blue", mfc="none", label=r"exact atoms ($\lambda\in\{0,1\}$)")
    for n in N[~pos]:
        ax[0].annotate("no atoms", (n + 0.28, 1.3e-3), ha="center", va="bottom", fontsize=7.5, color="tab:blue", rotation=90)
    ax[0].set_yscale("log")
    ax[0].set_ylim(1e-3, 1.5)
    ax[0].set_xlabel(r"$N'$")
    ax[0].set_ylabel(r"fraction with $\min(\lambda,1-\lambda)<10^{-3}$")
    ax[0].set_title("walls: arithmetic hard edge", fontsize=10)
    ax[0].axvspan(14.5, 16.5, color="0.92", zorder=0)
    ax[0].text(15.5, 1.05, r"$r>1$", ha="center", va="bottom", fontsize=8.5)
    ax[0].legend(frameon=False, loc="center left", bbox_to_anchor=(0.0, 0.45))
    axr = ax[0].twinx()
    axr.plot(N, [r["r"] for r in rows], color="tab:green", lw=0.9, alpha=0.7)
    axr.axhline(1.0, color="tab:green", lw=0.6, ls=":", alpha=0.7)
    axr.set_ylabel(r"$r=e_{F+2}/b_F$", color="tab:green")
    axr.tick_params(axis="y", colors="tab:green")
    axr.set_ylim(0, 2.2)
    NL = [n for n in N if n in het and "L_meas" in het[n]]
    ax[1].errorbar(NL, [het[n]["L_meas"] for n in NL], yerr=[het[n]["L_meas_se"] for n in NL], fmt="ko-", ms=5,
                   label=r"measured: $T_4^{\rm dec}-T_4^{\rm OU}(x_*)$")
    ax[1].plot(NL, [het[n]["L_pred"] for n in NL], "D--", color="tab:red", mfc="none", ms=5,
               label=rf"predicted: $T_2^{{\rm dec}}-{a.sym}_{{\rm first}}-T_4^{{\rm OU}}(x_*)$")
    ax[1].axvspan(14.5, 15.5, color="0.92", zorder=0)
    ax[1].set_ylim(0, 0.06)
    ax[1].set_xlim(7.5, 16.5)
    ax[1].set_xlabel(r"$N'$")
    ax[1].set_ylabel(r"localization $L$")
    ax[1].set_title("localization does not follow the walls", fontsize=10)
    ax[1].legend(frameon=False, loc="lower right")
    fig.tight_layout()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.savefig(a.out)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
