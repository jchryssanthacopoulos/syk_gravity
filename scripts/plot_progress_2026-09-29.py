#!/usr/bin/env python3
r"""
Figures for research/pdfs/progress_2026-09-29.pdf (data: results/data/parity_family_2026-09-29/).

  fig_parity_family      r/a and T_4/T_4^free vs parity-set size s (N' = 12, 13, 14; SYK and Haar)
  fig_size_transmission  (a) size spectrum w_k of the BPS projector, SYK vs Haar, with BLY's P_n at k = 2n (N'=14)
                         (b) relative transmission excess phi_k/a^2 - 1 vs size k (N' = 12, 13, 14; Haar at 14)
  fig_T4_decomposition   decoder T_4 - T_4^free split exactly into variance / two-model / remainder R
Light (< 0.5 GB).  Usage: .venv/bin/python scripts/plot_progress_2026-09-29.py [--figs DIR]
"""
import argparse
import csv
import json
import os
import sys

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from dssyk_n2 import hh_length_distribution  # noqa: E402
from nullmodels import free_compression_mu_moments  # noqa: E402
from size_decomposition import transmission_matrix, transmissions_from_weights  # noqa: E402

DATA = os.path.join(ROOT, "results/data/parity_family_2026-09-29")
BLUE, ORANGE, AQUA, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#8a8984"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
NCOL = {12: BLUE, 13: ORANGE, 14: AQUA}


def style(plt):
    plt.rcParams.update({"font.size": 10.5, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                         "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white"})


def fig_parity(plt, figs):
    rows = list(csv.DictReader(open(os.path.join(DATA, "summary.csv"))))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 3.9))
    for N in (12, 13, 14):
        for model, mk, col, ls in (("syk", "o", NCOL[N], "-"), ("haar", "s", GREY, ":")):
            R = sorted([r for r in rows if int(r["Np"]) == N and r["model"] == model], key=lambda r: int(r["s"]))
            s = [int(r["s"]) for r in R]
            ra = [float(r["r"]) / float(r["a"]) for r in R]
            t4 = [float(r["T4_over_free"]) for r in R]
            t4e = [float(r["T4_se"]) if r["T4_se"] not in ("nan", "") else 0 for r in R]
            lab = f"SYK N′={N}" if model == "syk" else ("Haar (all N′)" if N == 14 else None)
            ax1.plot(s, ra, marker=mk, ls=ls, color=col, ms=5.5, lw=1.4, label=lab, mec="white", mew=0.8)
            ax2.errorbar(s, t4, yerr=t4e, marker=mk, ls=ls, color=col, ms=5.5, lw=1.4, label=lab, mec="white",
                         mew=0.8, capsize=0)
    for ax, lab, title in ((ax1, r"$r/a$  (second moment / free)", "(a) variance → free as $x\\to0$"),
                           (ax2, r"$T_4/T_4^{\rm free}$  (fourth moment / free)", "(b) fourth moment does not")):
        ax.axhline(1, color=INK2, lw=0.9)
        ax.set_xlabel(r"parity-set size $s$   ($s=1$: the decoder;  $s=N'/2$: $x=0$)")
        ax.set_ylabel(lab)
        ax.set_title(title, fontsize=10.5)
        ax.grid(color=GRID, lw=0.7)
        ax.set_xticks(range(1, 8))
    ax1.legend(fontsize=8.5, frameon=False)
    ax2.set_ylim(0.95, 2.0)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(figs, f"fig_parity_family.{ext}"), dpi=170)


def fig_size(plt, figs):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 3.9))
    N = 14
    rows = [json.loads(l) for l in open(os.path.join(DATA, f"size_q3_N{N}.jsonl"))]
    wS = np.mean([r["w"] for r in rows if r["model"] == "syk"], axis=0)
    wH = np.mean([r["w"] for r in rows if r["model"] == "haar"], axis=0)
    k = np.arange(len(wS))
    ax1.bar(k - 0.2, wS, width=0.4, color=BLUE, label="SYK BPS projector (4 real.)")
    ax1.bar(k + 0.2, wH, width=0.4, color=GREY, label="Haar projector, same rank")
    Pn = hh_length_distribution(18 / N, 0.0)
    nmax = len(wS[0::2])
    ax1.plot(2 * np.arange(nmax), Pn[:nmax], "D", color=ORANGE, ms=7, mec="white", mew=1,
             label=r"BLY wormhole length $P_n$ at $k=2n$")
    ax1.set_xticks(k)
    ax1.set_xlabel(r"operator size $k$  (sector $V_k$ of End$(\Lambda^F)$)")
    ax1.set_ylabel(r"size weight $w_k=\|P^{(k)}\|^2/d$")
    ax1.set_title(r"(a) size spectrum of $P$, $N'=14$  ($w_0=a$ exactly)", fontsize=10.5)
    ax1.grid(axis="y", color=GRID, lw=0.7)
    ax1.legend(fontsize=8.5, frameon=False)
    for N, col in ((12, BLUE), (13, ORANGE), (14, AQUA)):
        P = N // 2
        A = transmission_matrix(N, P)
        R = [json.loads(l) for l in open(os.path.join(DATA, f"size_q3_N{N}.jsonl"))]
        for model, mk, ls, c in (("syk", "o", "-", col), ("haar", "s", ":", GREY)):
            RR = [r for r in R if r["model"] == model]
            if not RR:
                continue
            a, d = RR[0]["a"], RR[0]["d"]
            phi = np.mean([transmissions_from_weights(np.array(r["w"]), d, N, P, A) for r in RR], axis=0)
            kk = np.arange(1, len(phi))
            lab = f"SYK N′={N}" if model == "syk" else ("Haar" if N == 14 else None)
            ax2.plot(kk, phi[1:] / a ** 2 - 1, marker=mk, ls=ls, color=c, ms=5.5, lw=1.4, label=lab, mec="white",
                     mew=0.8)
    ax2.axhline(0, color=INK2, lw=0.9)
    ax2.set_xlabel(r"operator size $k$")
    ax2.set_ylabel(r"$\phi_k/a^2-1$  (transmission excess)")
    ax2.set_title(r"(b) transmission $X\mapsto E[PXP]$ on $V_k$: free $\Leftrightarrow$ flat", fontsize=10.5)
    ax2.grid(color=GRID, lw=0.7)
    ax2.legend(fontsize=8.5, frameon=False)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(figs, f"fig_size_transmission.{ext}"), dpi=170)


def fig_T4(plt, figs):
    Ns, lin, twm, rem, tot = [], [], [], [], []
    for N in (12, 13, 14):
        rows = [json.loads(l) for l in open(os.path.join(DATA, f"decoderT4_q3_N{N}.jsonl")) if json.loads(l)["model"] == "syk"]
        a = rows[0]["a"]
        b = 1 - (N // 2) / N
        fm = free_compression_mu_moments(a, b, 4)
        T2f, T4f = fm[1], fm[3]
        Q4f = T4f - a * a - 2 * a * (T2f - a)
        L, Tw, Rm, Tt = [], [], [], []
        for r in rows:
            S = r["c0"] - a * a
            for e in r["probes"]:
                if e["kind"] == "site":
                    L.append(2 * a * (e["T2"] - T2f)); Tw.append(S - Q4f); Rm.append(e["R_exact"]); Tt.append(e["T4"] - T4f)
        Ns.append(N); lin.append(np.mean(L)); twm.append(np.mean(Tw)); rem.append(np.mean(Rm)); tot.append(np.mean(Tt))
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    x = np.arange(len(Ns))
    ax.bar(x, lin, width=0.55, color=BLUE, label=r"one nontrivial $P'$: $2a(T_2-T_2^{\rm free})$ (exact from variance)")
    ax.bar(x, twm, width=0.55, bottom=lin, color=AQUA, label=r"two-model part $S-Q_4^{\rm free}$ (exact from $w_k$)")
    ax.bar(x, rem, width=0.55, bottom=np.array(lin) + np.array(twm), color=ORANGE,
           label=r"remainder $R$ (four-point, open)")
    for i in range(len(Ns)):
        ax.text(x[i], tot[i] + 0.006, f"{rem[i] / tot[i] * 100:.0f}% R", ha="center", fontsize=9, color=INK2)
    ax.set_xticks(x)
    ax.set_xticklabels([f"N′={n}" for n in Ns])
    ax.set_ylabel(r"$T_4-T_4^{\rm free}$  (decoder)")
    ax.set_title("Exact decomposition of the decoder's fourth-moment excess", fontsize=10.5)
    ax.grid(axis="y", color=GRID, lw=0.7)
    ax.set_ylim(0, max(tot) * 1.42)
    ax.legend(fontsize=8.5, frameon=False, loc="upper left")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(figs, f"fig_T4_decomposition.{ext}"), dpi=170)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figs", default=os.path.join(ROOT, "results/figures/progress_2026-09-29"))
    a = ap.parse_args()
    os.makedirs(a.figs, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    style(plt)
    fig_parity(plt, a.figs)
    fig_size(plt, a.figs)
    fig_T4(plt, a.figs)
    print("wrote", a.figs)


if __name__ == "__main__":
    main()
