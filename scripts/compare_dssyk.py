#!/usr/bin/env python3
r"""
Compare exact q=3 SYK decoder data with the finite-lambda super-chord predictions of Boruch-Lin-Yan
(src/dssyk_n2.py; docs/derivations.md D2) and with the super-Schwarzian (LMRS; src/lmrs_predictions.py).
No free parameters: lambda = 2 p^2 / N' = 18/N', Delta = 1/3, j = (P - N'/2)/3.

Reads   results/data/chord_invariants_2026-09-28/summary.csv   (r, a; N' <= 15)
        results/data/chord_invariants_2026-09-28/raw_q3_N16.jsonl (N' = 16, if present)
        results/data/lmrs_variance_2026-09-28/jsector.jsonl       (j = +-1/3 sectors)
Writes  results/data/dssyk_2026-09-28/comparison.csv
        results/figures/dssyk_2026-09-28/fig_r_a_vs_chords.{png,pdf}
Light (< 0.2 GB).
"""
import argparse
import csv
import json
import os
import sys

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from dssyk_n2 import bps_fraction_chord, bps_fraction_index, r_chord  # noqa: E402
from lmrs_predictions import r_pred  # noqa: E402

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=os.path.join(ROOT, "results/data/dssyk_2026-09-28"))
    ap.add_argument("--figs", default=os.path.join(ROOT, "results/figures/dssyk_2026-09-28"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    os.makedirs(a.figs, exist_ok=True)
    data = os.path.join(ROOT, "results/data/chord_invariants_2026-09-28")

    rows = {int(r["Np"]): (float(r["r"]), float(r["r_se"])) for r in csv.DictReader(open(os.path.join(data, "summary.csv")))
            if r["q"] == "3" and r["model"] == "syk"}
    f16 = os.path.join(data, "raw_q3_N16.jsonl")
    if os.path.exists(f16):                       # same reduction as analyze_chord_invariants: mean r over modes
        rr = [np.mean([s["rho"][1] for s in json.loads(l)["single"]]) for l in open(f16)
              if json.loads(l)["model"] == "syk"]
        if rr:
            rows[16] = (float(np.mean(rr)), float(np.std(rr, ddof=1) / np.sqrt(len(rr))) if len(rr) > 1 else np.nan)
    js = [json.loads(l) for l in open(os.path.join(ROOT, "results/data/lmrs_variance_2026-09-28/jsector.jsonl"))]

    out = []
    for N in sorted(rows):
        P = N // 2
        r, se = rows[N]
        j = (P - N / 2) / 3
        out.append(dict(Np=N, P=P, j=j, lam=18 / N, r=r, r_se=se, r_chord=r_chord(N), r_lmrs=r_pred(N),
                        a_exact=bps_fraction_index(N, P), a_chord=bps_fraction_chord(18 / N, j)))
    for N in sorted({x["Np"] for x in js}):
        P = N // 2 - 1
        b = 1 - P / N
        v = [x["sectors"][str(P)]["var"] / (b * (1 - b)) for x in js if x["Np"] == N]
        out.append(dict(Np=N, P=P, j=(P - N / 2) / 3, lam=18 / N, r=float(np.mean(v)),
                        r_se=float(np.std(v, ddof=1) / np.sqrt(len(v))), r_chord=r_chord(N, P), r_lmrs=r_pred(N, P),
                        a_exact=bps_fraction_index(N, P), a_chord=bps_fraction_chord(18 / N, (P - N / 2) / 3)))
    with open(os.path.join(a.out, "comparison.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print("N'   j       r (SYK)          r_chord  r_LMRS  | meas/chord  meas/LMRS |  a_exact  a_chord")
    for o in out:
        print(f"{o['Np']:3d} {o['j']:+.3f}  {o['r']:.4f}±{o['r_se']:.4f}  {o['r_chord']:.4f}  {o['r_lmrs']:.4f}  |"
              f"   {o['r'] / o['r_chord']:.3f}      {o['r'] / o['r_lmrs']:.3f}   |  {o['a_exact']:.4f}   {o['a_chord']:.4f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                         "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "savefig.facecolor": "#fcfcfb"})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.4))
    main_ = [o for o in out if abs(o["j"]) < 0.3]
    side = [o for o in out if abs(o["j"]) > 0.3]
    for par, mk, lab in ((0, "o", "even N′ (j = 0)"), (1, "s", "odd N′ (j = ∓1/6)")):
        sel = [o for o in main_ if o["Np"] % 2 == par]
        Ns = np.arange(6 + par, 21, 2)
        ax1.plot(Ns, [r_pred(n) for n in Ns], ls=":", marker=mk, ms=4, mfc="none", color=ORANGE, lw=1.4,
                 label=f"super-Schwarzian (LMRS) — {lab}")
        ax1.plot(Ns, [r_chord(n) for n in Ns], ls="-", marker=mk, ms=4, mfc="none", color=INK2, lw=1.4,
                 label=f"super-chord, finite λ (BLY) — {lab}")
        ax1.errorbar([o["Np"] for o in sel], [o["r"] for o in sel], yerr=[o["r_se"] for o in sel], fmt=mk,
                     color=BLUE, ms=7, mec="#fcfcfb", mew=1.2, zorder=5, label=f"SYK exact — {lab}")
    Ns = np.arange(8, 21, 2)
    ax1.plot(Ns, [r_chord(n, n // 2 - 1) for n in Ns], ls="-", marker="D", ms=4, mfc="none", color=AQUA, lw=1.4,
             label="super-chord, j = ±1/3 sector")
    ax1.errorbar([o["Np"] for o in side], [o["r"] for o in side], yerr=[o["r_se"] for o in side], fmt="D",
                 color=AQUA, ms=7, mec="#fcfcfb", mew=1.2, zorder=5, label="SYK exact, j = ±1/3 sector")
    ax1.set_xlabel("N′")
    ax1.set_ylabel("r = Var(λ) / [b(1−b)]")
    ax1.set_title("(a) decoder variance ratio: no free parameters", fontsize=11)
    ax1.set_ylim(0.3, 1.4)
    ax1.grid(color=GRID, lw=0.7)
    ax1.legend(fontsize=7, frameon=False, loc="upper right", ncol=1)

    x = [o["Np"] for o in main_]
    ax2.plot(x, [o["r"] / o["r_chord"] for o in main_], "o", color=BLUE, ms=7, mec="#fcfcfb", mew=1.2,
             label="r: SYK / super-chord (j = 0, ∓1/6)")
    ax2.plot([o["Np"] for o in side], [o["r"] / o["r_chord"] for o in side], "D", color=AQUA, ms=7, mec="#fcfcfb",
             mew=1.2, label="r: SYK / super-chord (j = ±1/3)")
    ax2.plot(x, [o["r"] / o["r_lmrs"] for o in main_], "o", color=ORANGE, ms=7, mec="#fcfcfb", mew=1.2,
             label="r: SYK / super-Schwarzian (j = 0, ∓1/6)")
    ax2.plot(x, [o["a_exact"] / o["a_chord"] for o in main_], "^", color=INK2, ms=6, mfc="none",
             label="BPS fraction a: exact / super-chord")
    ax2.axhline(1, color=INK2, lw=1)
    ax2.set_xlabel("N′")
    ax2.set_ylabel("exact / prediction")
    ax2.set_title("(b) ratios to the predictions", fontsize=11)
    ax2.set_ylim(0.6, 1.12)
    ax2.grid(color=GRID, lw=0.7)
    ax2.legend(fontsize=8, frameon=False, loc="lower right")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(a.figs, f"fig_r_a_vs_chords.{ext}"), dpi=160)
    print("wrote", a.figs)


if __name__ == "__main__":
    main()
