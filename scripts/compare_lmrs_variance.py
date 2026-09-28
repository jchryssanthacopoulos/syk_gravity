#!/usr/bin/env python3
r"""
Compare the measured decoder variance ratio r = Var(lambda)/[b(1-b)] (chord-invariant campaign, q = 3) with the
parameter-free LMRS zero-energy prediction (src/lmrs_predictions.py; derivation docs/derivations.md D1).

Reads  results/data/chord_invariants_2026-09-28/summary.csv  (written by scripts/analyze_chord_invariants.py).
Writes results/data/lmrs_variance_2026-09-28/comparison.csv, fit.json
       results/figures/lmrs_variance_2026-09-28/fig_r_vs_lmrs.{png,pdf}
Fits (weighted by the s.e. of r, N' >= --nmin): ratio = 1 - c/N'  and  ratio = A - B/N'.
Light computation (< 0.2 GB); no memory watchdog needed.

Usage: .venv/bin/python scripts/compare_lmrs_variance.py [--nmin 9] [--alpha-s 0.00842]
"""
import argparse
import csv
import json
import os
import sys

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from lmrs_predictions import ALPHA_S_Q3, r_charge, r_pred  # noqa: E402

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"     # same validated slots as analyze_chord_invariants.py
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary", default=os.path.join(ROOT, "results/data/chord_invariants_2026-09-28/summary.csv"))
    ap.add_argument("--nmin", type=int, default=9)
    ap.add_argument("--alpha-s", type=float, default=ALPHA_S_Q3)
    ap.add_argument("--out", default=os.path.join(ROOT, "results/data/lmrs_variance_2026-09-28"))
    ap.add_argument("--figs", default=os.path.join(ROOT, "results/figures/lmrs_variance_2026-09-28"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    os.makedirs(a.figs, exist_ok=True)

    rows = [r for r in csv.DictReader(open(a.summary)) if r["q"] == "3" and r["model"] == "syk"]
    N = np.array([int(r["Np"]) for r in rows])
    R = np.array([float(r["r"]) for r in rows])
    SE = np.array([float(r["r_se"]) for r in rows])
    A_ = np.array([float(r["a"]) for r in rows])
    n_real = [int(r["n"]) for r in rows]
    PR = np.array([r_pred(n, alpha_S=a.alpha_s) for n in N])
    rat, rse = R / PR, SE / PR

    m = N >= a.nmin
    w = 1 / rse[m] ** 2
    c = float(np.sum(w * (1 - rat[m]) / N[m]) / np.sum(w / N[m] ** 2))
    Xd = np.vstack([np.ones(m.sum()), -1 / N[m]]).T * np.sqrt(w)[:, None]
    (A2, B2), *_ = np.linalg.lstsq(Xd, rat[m] * np.sqrt(w), rcond=None)
    fit = dict(alpha_S=a.alpha_s, nmin=a.nmin, N_fit=N[m].tolist(),
               one_param=dict(c=c, max_abs_resid=float(np.max(np.abs(rat[m] - (1 - c / N[m]))))),
               two_param=dict(A=float(A2), B=float(B2),
                              max_abs_resid=float(np.max(np.abs(rat[m] - (A2 - B2 / N[m]))))),
               predictions={str(n): dict(r_LMRS=r_pred(n, alpha_S=a.alpha_s),
                                         r_1param=r_pred(n, alpha_S=a.alpha_s) * (1 - c / n),
                                         r_2param=r_pred(n, alpha_S=a.alpha_s) * (A2 - B2 / n))
                            for n in (16, 17, 18, 20, 24, 32)})
    # Alternative models (is the LMRS N'^{-2/3} shape actually preferred?  see docs/derivations.md D1)
    from scipy.optimize import least_squares
    Nm, Rm, Sm, Pm, Am = N[m], R[m], SE[m], PR[m], A_[m]
    odd = (Nm % 2).astype(float)
    jfac = np.where(odd > 0, np.array([r_pred(n, alpha_S=a.alpha_s) / r_pred(n - 1, alpha_S=a.alpha_s)
                                       * ((n - 1) / n) ** (-2 / 3) for n in Nm]), 1.0)   # parity part of LMRS
    models = {
        "LMRS (0 params)": (0, lambda p: Pm, []),
        "LMRS x A (= alpha_S free)": (1, lambda p: p[0] * Pm, [0.8]),
        "LMRS x (1 - c/N)": (1, lambda p: Pm * (1 - p[0] / Nm), [2.0]),
        "K N^-p (no parity)": (2, lambda p: p[0] * Nm ** (-p[1]), [2.0, 0.5]),
        "K N^-p x LMRS parity factor": (2, lambda p: p[0] * Nm ** (-p[1]) * jfac, [2.0, 0.5]),
        "K N^-p x (1 - e odd)": (3, lambda p: p[0] * Nm ** (-p[1]) * (1 - p[2] * odd), [2.0, 0.5, 0.05]),
        "Haar a x (1 + c/N)": (1, lambda p: Am * (1 + p[0] / Nm), [1.0]),
    }
    alt = {}
    for name, (k, f, x0) in models.items():
        x = least_squares(lambda p: (f(p) - Rm) / Sm, x0).x if k else np.array([])
        res = f(x) - Rm
        alt[name] = dict(n_params=k, params=x.tolist(), rms=float(np.sqrt(np.mean(res ** 2))),
                         chi2_per_dof=float(np.sum((res / Sm) ** 2) / (len(Nm) - k)))
        print(f"  {name:32s} k={k} params={np.round(x, 3)} rms={alt[name]['rms']:.4f} "
              f"chi2/dof={alt[name]['chi2_per_dof']:.1f}")
    fit["alternative_models"] = alt
    with open(os.path.join(a.out, "fit.json"), "w") as f:
        json.dump(fit, f, indent=1)
    with open(os.path.join(a.out, "comparison.csv"), "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["Np", "j", "n_real", "a", "r", "r_se", "r_LMRS", "ratio", "ratio_se"])
        for i, n in enumerate(N):
            wr.writerow([n, r_charge(n, n // 2, 3), n_real[i], A_[i], R[i], SE[i], PR[i], rat[i], rse[i]])

    print(f"alpha_S = {a.alpha_s}; fit N' >= {a.nmin}")
    print(f"ratio = 1 - c/N':  c = {c:.3f}  (max |resid| {fit['one_param']['max_abs_resid']:.4f})")
    print(f"ratio = A - B/N':  A = {A2:.3f}, B = {B2:.3f}  (max |resid| {fit['two_param']['max_abs_resid']:.4f})")
    for n, p in fit["predictions"].items():
        print(f"  N'={n}: r_LMRS={p['r_LMRS']:.4f}  x(1-c/N)={p['r_1param']:.4f}  x(A-B/N)={p['r_2param']:.4f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                         "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
                         "savefig.facecolor": "#fcfcfb"})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.4))
    ev, od = N % 2 == 0, N % 2 == 1
    Nc = np.linspace(7, 20, 300)
    for par, mk, lab in ((0, "o", "even N′ (j = 0)"), (1, "s", "odd N′ (j = −1/6)")):
        sel = N % 2 == par
        Ns = np.arange(8 if par == 0 else 7, 21, 2)
        ax1.plot(Ns, [r_pred(n, alpha_S=a.alpha_s) for n in Ns], ls=":", color=ORANGE, lw=1.5, marker=mk,
                 ms=4, mfc="none", label=f"LMRS, no free parameters — {lab}")
        ax1.plot(Ns, [r_pred(n, alpha_S=a.alpha_s) * (1 - c / n) for n in Ns], ls="-", color=ORANGE, lw=2,
                 alpha=0.55)
        ax1.errorbar(N[sel], R[sel], yerr=SE[sel], fmt=mk, color=BLUE, ms=7, mec="#fcfcfb", mew=1.2,
                     label=f"SYK exact — {lab}", zorder=5)
    ax1.plot(N, A_, "^", color=AQUA, ms=7, mec="#fcfcfb", mew=1.2, label="Haar / Wachter: r = a")
    ax1.plot([], [], color=ORANGE, lw=2, alpha=0.55, label=f"LMRS × (1 − {c:.2f}/N′), c fitted")
    ax1.axhline(1, color=GRID, lw=1)
    ax1.set_xlabel("N′ (fermions in the enlarged theory)")
    ax1.set_ylabel("r = Var(λ) / [b(1−b)]")
    ax1.set_title("(a) decoder variance ratio vs super-Schwarzian", fontsize=11)
    ax1.set_ylim(0.2, 1.2)
    ax1.set_xlim(6.5, 20.5)
    ax1.grid(color=GRID, lw=0.7)
    ax1.legend(fontsize=8, frameon=False, loc="lower left")

    x = 1 / N
    ax2.errorbar(x[ev], rat[ev], yerr=rse[ev], fmt="o", color=BLUE, ms=7, mec="#fcfcfb", mew=1.2,
                 label="even N′ (j = 0)")
    ax2.errorbar(x[od], rat[od], yerr=rse[od], fmt="s", color=BLUE, ms=7, mec="#fcfcfb", mew=1.2,
                 label="odd N′ (j = −1/6)")
    xc = np.linspace(0, 1 / 5, 100)
    ax2.plot(xc, 1 - c * xc, color=ORANGE, lw=2, label=f"1 − {c:.2f}/N′ (c fitted, N′ ≥ {a.nmin})")
    ax2.plot(xc, A2 - B2 * xc, color=INK2, lw=1.2, ls="--", label=f"{A2:.3f} − {B2:.2f}/N′")
    ax2.axhline(1, color=GRID, lw=1)
    ax2.set_xlim(0, 0.21)
    ax2.set_ylim(0.55, 1.05)
    ax2.set_xlabel("1/N′")
    ax2.set_ylabel("r_SYK / r_LMRS")
    ax2.set_title("(b) ratio to the parameter-free prediction", fontsize=11)
    ax2.grid(color=GRID, lw=0.7)
    ax2.legend(fontsize=8.5, frameon=False, loc="lower left")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(a.figs, f"fig_r_vs_lmrs.{ext}"), dpi=160)
    print("wrote", a.figs)


if __name__ == "__main__":
    main()
