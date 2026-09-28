#!/usr/bin/env python3
r"""
Figure: "deviation from freeness = wormholes of nonzero length" (docs/derivations.md D2, Implications).

(a) Chord-number distribution P_n of the supersymmetric Hartle-Hawking state (BLY) at N' = 14 (lambda = 18/14,
    j = 0), and its contribution P_n q^{2 Delta n} to r (Delta = 1/3).  The n = 0 bar is the BPS fraction a, i.e.
    the double-scaled BPS fraction a_lambda = D(j), which corresponds to the Haar / free-compression value of r
    (equality of values within the double-scaled theory; exact a is 0-10% smaller at p = 3).
(b) Rigidity fraction eta_r = (r - a)/(1 - a) = <q^{2 Delta n}>_{n>=1} predicted by the chord theory vs N',
    with the measured SYK values (chord-invariant campaign).
No free parameters.  Light (< 0.2 GB).

Usage: .venv/bin/python scripts/plot_wormhole_freeness.py [--figs DIR]
"""
import argparse
import csv
import os
import sys
from math import exp

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from dssyk_n2 import bps_fraction_chord, hh_length_distribution, two_point_chord  # noqa: E402

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figs", default=os.path.join(ROOT, "results/figures/dssyk_2026-09-28"))
    a = ap.parse_args()
    os.makedirs(a.figs, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                         "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white"})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.2))

    N, lam = 14, 18 / 14
    P = hh_length_distribution(lam, 0.0)
    n = np.arange(len(P))
    contrib = P * exp(-lam) ** (2 * n / 3)
    k = 12
    ax1.bar(n[:k] - 0.2, P[:k], width=0.4, color=GRID, edgecolor=INK2, lw=0.6, label=r"$P_n$ (wormhole length)")
    ax1.bar(n[:k] + 0.2, contrib[:k], width=0.4, color=BLUE, label=r"$P_n\,q^{2\Delta n}$ (contribution to $r$)")
    ax1.bar([0.2], [contrib[0]], width=0.4, color=ORANGE, label=r"$n=0$: $a_\lambda=P_0$ ($\leftrightarrow$ Haar/free value)")
    r, aa = two_point_chord(lam, 0.0, 1 / 3), bps_fraction_chord(lam, 0.0)
    ax1.text(5.3, 0.13, rf"$\sum_n P_n = 1$" "\n" rf"$r_\lambda=\sum_n P_n q^{{2\Delta n}} = {r:.3f}$" "\n"
             rf"$a_\lambda = P_0 = {aa:.3f}$" "\n" rf"$r_\lambda-a_\lambda = \sum_{{n\geq1}} = {r - aa:.3f}$", fontsize=10, color=INK)
    ax1.set_xlabel(r"chord number $n$ (wormhole length $\ell = 2\lambda n$)")
    ax1.set_ylabel("probability / contribution")
    ax1.set_title(rf"(a) SUSY wormhole at $N'={N}$ ($\lambda={lam:.2f}$, $j=0$, $\Delta=1/3$)", fontsize=11)
    ax1.set_xticks(range(0, k))
    ax1.legend(fontsize=8.5, frameon=False, loc="upper right")
    ax1.grid(axis="y", color=GRID, lw=0.7)

    Ns = np.arange(6, 101, 2)
    eta = [(two_point_chord(18 / m, 0.0, 1 / 3) - bps_fraction_chord(18 / m, 0.0)) /
           (1 - bps_fraction_chord(18 / m, 0.0)) for m in Ns]
    ax2.plot(Ns, eta, color=INK2, lw=2, label=r"super-chord, $j=0$: $(r_\lambda-a_\lambda)/(1-a_\lambda)$, no free parameters")
    rows = [r_ for r_ in csv.DictReader(open(os.path.join(ROOT, "results/data/chord_invariants_2026-09-28/summary.csv")))
            if r_["q"] == "3" and r_["model"] == "syk" and int(r_["Np"]) % 2 == 0]
    xm = [int(r_["Np"]) for r_ in rows]
    ym = [(float(r_["r"]) - float(r_["a"])) / (1 - float(r_["a"])) for r_ in rows]
    ax2.plot(xm, ym, "o", color=BLUE, ms=7, mec="white", mew=1.2, label=r"SYK exact, even $N'$: $(r-a)/(1-a)$ with exact $a$")
    ax2.axvline(22, color=GRID, lw=1, ls="--")
    ax2.text(23, 0.40, "predicted maximum\n" r"near $N'\approx 22$", fontsize=9, color=INK2)
    ax2.set_xscale("log")
    ax2.set_xticks([6, 10, 20, 50, 100])
    ax2.set_xticklabels(["6", "10", "20", "50", "100"])
    ax2.set_xlabel(r"$N'$")
    ax2.set_ylabel(r"$\eta_r=(r-a)/(1-a)=\langle q^{2\Delta n}\rangle_{n\geq1}$")
    ax2.set_title(r"(b) rigidity fraction vs $N'$", fontsize=11)
    ax2.set_ylim(0, 0.45)
    ax2.grid(color=GRID, lw=0.7)
    ax2.legend(fontsize=8.5, frameon=False, loc="lower center")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(a.figs, f"fig_wormhole_freeness.{ext}"), dpi=160)
    print("wrote", a.figs)


if __name__ == "__main__":
    main()
