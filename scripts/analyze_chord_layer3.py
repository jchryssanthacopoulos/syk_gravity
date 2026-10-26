#!/usr/bin/env python3
r"""
Analysis of results/data/chord_layer3_2026-10-26 (scripts/run_chord_layer3.py; docs/derivations.md D11 §4).

Table: seed means per N' of T2, H and the low-tail fractions P(lambda < x) for each R model:
  exact, haar, flat, cR_Ew, chord, diag, zero (= first order).
The exact rank fraction of R is compared with exp(-10 lambda/9).
Figure results/figures/chord_layer3_2026-10-26/fig_layer3.pdf:
  (a) low tail: cumulative distribution of lambda at N' = 12, 14, 16 for exact, first order, Haar, and the full chord
      model;
  (b) T2 and P(lambda < 0.1) against N' for each model;
  (c) rank fraction of R: exact against exp(-10 lambda/9).
Usage: .venv/bin/python scripts/analyze_chord_layer3.py
"""
import glob
import json
import os
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DATA = os.path.join(ROOT, "results/data/chord_layer3_2026-10-26")
FIGS = os.path.join(ROOT, "results/figures/chord_layer3_2026-10-26")
MODELS = ("exact", "haar", "flat", "cR_Ew", "chord", "diag", "zero")
plt.rcParams.update({"font.size": 9, "font.family": "serif", "mathtext.fontset": "cm"})
COL = {"exact": "k", "haar": "#1b9e77", "flat": "#66a61e", "cR_Ew": "#7570b3", "chord": "#d95f02", "diag": "#e7298a",
       "zero": "0.55"}
LAB = {"exact": "exact", "haar": "R Haar-rotated", "flat": "flat spectrum, Haar", "cR_Ew": r"$c_R/E_w$, Haar",
       "chord": "full chord model", "diag": "diag R (mean field)", "zero": "first order (R = 0)"}


def load():
    recs = []
    for fn in sorted(glob.glob(os.path.join(DATA, "layer3*.jsonl"))):
        recs += [json.loads(x) for x in open(fn)]
    return recs


def pooled(Np, key):
    xs = [np.load(fn)[key] for fn in sorted(glob.glob(os.path.join(DATA, f"lam_N{Np}_s*.npz")))]
    return np.sort(np.concatenate(xs)) if xs else np.zeros(0)


def main():
    os.makedirs(FIGS, exist_ok=True)
    recs = load()
    g = defaultdict(list)
    for r in recs:
        g[r["Np"]].append(r)
    rows = []
    print(" N' seeds  r_R/n (exp(-10lam/9)) | model: T2  H  P(<.01) P(<.03) P(<.1) min")
    for Np in sorted(g):
        rs = g[Np]
        row = dict(Np=Np, seeds=len(rs), rank_frac=float(np.mean([r["rank_frac"] for r in rs])),
                   rank_frac_chord=rs[0]["rank_frac_chord"], cR=float(np.mean([r["cR"] for r in rs])),
                   cR_chord=rs[0]["cR_chord"], rho_mean=float(np.mean([r["rho_mean"] for r in rs])),
                   rho_chord_mean=float(np.mean([r["rho_chord_mean"] for r in rs])))
        print(f"{Np:3d}  {len(rs)}   {row['rank_frac']:.3f} ({row['rank_frac_chord']:.3f})  c_R {row['cR']:.1f} (chord {row['cR_chord']:.1f})"
              f"  rho mean {row['rho_mean']:.3f} (chord {row['rho_chord_mean']:.3f})")
        for m in MODELS:
            vals = {k: float(np.mean([r["models"][m][k] for r in rs])) for k in ("T2", "H", "P_lt_0.01", "P_lt_0.03", "P_lt_0.1")}
            vals["lam_min"] = float(np.min([r["models"][m]["lam_min"] for r in rs]))
            row[m] = vals
            print(f"      {m:6s}: {vals['T2']:.4f} {vals['H']:.4f} {vals['P_lt_0.01']:.3f} {vals['P_lt_0.03']:.3f}"
                  f" {vals['P_lt_0.1']:.3f} {vals['lam_min']:.2e}")
        rows.append(row)
    json.dump(rows, open(os.path.join(DATA, "summary.json"), "w"), indent=1)

    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6))
    ax = axs[0]
    labelled = set()
    for Np, ls in ((12, ":"), (14, "--"), (16, "-")):
        for m in ("exact", "zero", "haar", "chord"):
            x = pooled(Np, m)
            if not len(x):
                continue
            cc = np.arange(1, len(x) + 1) / len(x)
            ax.semilogx(np.maximum(x, 1e-7), cc, ls=ls, color=COL[m], lw=1.1 if m != "exact" else 1.4,
                        label=LAB[m] if m not in labelled else None)
            labelled.add(m)
    ax.set_xlim(1e-6, 1)
    ax.set_ylim(0, 0.2)
    ax.set_xlabel(r"decoder eigenvalue $\lambda$ (resonant $\ell'=-2$)")
    ax.set_ylabel("cumulative fraction")
    ax.legend(fontsize=5.8, frameon=False, loc="upper left", title=r"$N'=12$ (:), 14 (- -), 16 (—)", title_fontsize=5.8)
    ax.set_title("(a) low tail", fontsize=8.5, loc="left")
    ax = axs[1]
    Ns = [r["Np"] for r in rows]
    for m in MODELS:
        ax.plot(Ns, [r[m]["T2"] for r in rows], "o-", color=COL[m], ms=3.5, lw=1.0, label=LAB[m])
    ax.set_xlabel(r"$N'$")
    ax.set_ylabel(r"$T_2$")
    ax.legend(fontsize=5.4, frameon=False, loc="lower left")
    ax.set_title(r"(b) $T_2$ by model", fontsize=8.5, loc="left")
    ax = axs[2]
    ax.plot(Ns, [r["rank_frac"] for r in rows], "ko", ms=4, label="exact rank R / n")
    xx = np.linspace(min(Ns) - 1, max(Ns) + 1, 50)
    ax.plot(xx, np.exp(-(10 / 9) * 18 / (xx - 1)), "-", color=COL["chord"], lw=1.1, label=r"$e^{-10\lambda/9}$, $\lambda=18/N$")
    ax.set_xlabel(r"$N'$")
    ax.set_ylabel(r"rank fraction of $R$")
    ax.legend(fontsize=6.5, frameon=False)
    ax.set_title("(c) back-reaction channel", fontsize=8.5, loc="left")
    fig.tight_layout(pad=0.4, w_pad=0.6)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"fig_layer3.{ext}"), dpi=200)
    print("wrote", os.path.join(DATA, "summary.json"), "and", os.path.join(FIGS, "fig_layer3.pdf"))


if __name__ == "__main__":
    main()
