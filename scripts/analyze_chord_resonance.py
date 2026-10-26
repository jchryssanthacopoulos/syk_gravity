#!/usr/bin/env python3
r"""
Analysis of results/data/chord_resonance_2026-10-26 (scripts/run_chord_resonance.py): layer 1 of the resonant chord
problem, i.e. the first-order dressing spectrum of the resonant family from the chord rules against finite N.

Error budget (as in D7 sec. 4), per sector and family, from seed means of the first-order decoder statistics
W = <y(1-y)> and T2 = <y>, with y = ((1 - sigma^2)/(1 + sigma^2))^2:
  planar error     wishart (exact inner S, factorized BPS faces) vs sampled (exact first order over the eta ensemble)
  continuum error  chord (D7 jump law at the sector's charges) vs wishart (exact S)
  total            chord vs sampled
Tails (pooled spectra from the .npz files):
  P(sigma^2 > K) sqrt(K) at K = 10, 100, 1000 and Hill exponents, for sampled, wishart and chord.

Figures in results/figures/chord_resonance_2026-10-26/:
  fig_layer1.pdf
  (a) CCDF of sigma^2 in the resonant sector at N' = 12, 14, 16: sampled, wishart, chord;
  (b) first-order decoder density lambda_1 = 1/(1 + sigma^2): sampled vs chord, N' = 14, 16, with the exact decoder
      for comparison;
  (c) relative error of W and T2 (chord vs sampled) against N', resonant vs central.

Usage: .venv/bin/python scripts/analyze_chord_resonance.py
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
DATA = os.path.join(ROOT, "results/data/chord_resonance_2026-10-26")
WIN = os.path.join(ROOT, "results/data/uplift_window_2026-10-26")
FIGS = os.path.join(ROOT, "results/figures/chord_resonance_2026-10-26")
plt.rcParams.update({"font.size": 9, "font.family": "serif", "mathtext.fontset": "cm"})
LAYERS = ("own", "sampled", "kron", "wishart", "chord")


def load():
    recs = []
    for fn in sorted(glob.glob(os.path.join(DATA, "chord_resonance*.jsonl"))):
        recs += [json.loads(x) for x in open(fn)]
    return recs


def pooled(Np, F, fam, key):
    xs = []
    for fn in sorted(glob.glob(os.path.join(DATA, f"spec_N{Np}_F{F}_{fam}_s*.npz"))):
        z = np.load(fn)
        if key in z:
            xs.append(z[key])
    return np.concatenate(xs) if xs else np.zeros(0)


def wy(s):
    s = np.clip(np.asarray(s, float), 0, None)
    y = ((1 - s) / (1 + s)) ** 2
    return y * (1 - y), y


def exact_decoder(Np, lp):
    """Exact (full uplift) h-family decoder eigenvalues of the window runs, pooled over seeds."""
    lam = []
    for fn in [os.path.join(WIN, "sectors_p3.jsonl")] + glob.glob(os.path.join(WIN, "sectors_p3_*.jsonl")):
        for line in open(fn):
            r = json.loads(line)
            if r["Np"] == Np and r["lp"] == lp and r["d_h"] and r["d_g"] <= 1:
                lam += r["lam"]
    return np.array(lam)


def table(recs):
    g = defaultdict(list)
    for r in recs:
        g[(r["Np"], r["lp"], r["fam"])].append(r)
    rows = []
    for (Np, lp, fam), rs in sorted(g.items()):
        m = lambda lay, k: float(np.mean([r[lay][k] for r in rs]))  # noqa: E731
        row = dict(Np=Np, lp=lp, fam=fam, seeds=len(rs), b=rs[0]["b"], e=rs[0]["e"],
                   Emin_rel=float(np.mean([r["Emin_rel"] for r in rs])))
        for lay in LAYERS:
            row[f"W_{lay}"], row[f"T2_{lay}"] = m(lay, "W"), m(lay, "T2")
        row["W_planar"] = row["W_wishart"] / row["W_sampled"] - 1
        row["W_continuum"] = row["W_chord"] / row["W_wishart"] - 1
        row["W_total"] = row["W_chord"] / row["W_sampled"] - 1
        row["T2_planar"] = row["T2_wishart"] / row["T2_sampled"] - 1
        row["T2_continuum"] = row["T2_chord"] / row["T2_wishart"] - 1
        row["T2_total"] = row["T2_chord"] / row["T2_sampled"] - 1
        F = rs[0]["F"]
        for lay in ("sampled", "wishart", "chord"):
            x = pooled(Np, F, fam, lay)
            if len(x):
                row[f"Psqrt100_{lay}"] = float(np.mean(x > 100) * 10)
                row[f"Psqrt1000_{lay}"] = float(np.mean(x > 1000) * np.sqrt(1000))
                t = x[x > 100]
                row[f"hill100_{lay}"] = float(len(t) / np.sum(np.log(t / 100))) if len(t) > 2 else float("nan")
                row[f"hill100_n_{lay}"] = int(len(t))
        rows.append(row)
    return rows


def figure(rows):
    os.makedirs(FIGS, exist_ok=True)
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6))
    ax = axs[0]
    cols = {12: "#1b9e77", 14: "#7570b3", 16: "#d95f02"}
    for Np in (12, 14, 16):
        F = (Np - 2) // 2
        for lay, ls in (("sampled", "-"), ("wishart", ":"), ("chord", "--")):
            x = np.sort(pooled(Np, F, "h", lay))
            if not len(x):
                continue
            cc = 1 - np.arange(len(x)) / len(x)
            k = x > 1e-8
            ax.loglog(x[k], cc[k], ls=ls, color=cols[Np], lw=1.0, label=f"N'={Np}" if lay == "sampled" else None)
    ax.plot([], [], "k-", lw=1, label="exact (η ensemble)")
    ax.plot([], [], "k:", lw=1, label="exact S, planar")
    ax.plot([], [], "k--", lw=1, label="chord")
    ax.set_xlim(1e-2, 1e7)
    ax.set_ylim(1e-4, 1.2)
    ax.set_xlabel(r"first-order dressing cost $K=\sigma^2$")
    ax.set_ylabel(r"$P(\sigma^2>K)$")
    ax.legend(fontsize=5.8, frameon=False, loc="lower left")
    ax.set_title(r"(a) resonant tail, $\ell'=-2$", fontsize=8.5, loc="left")
    ax = axs[1]
    bins = np.linspace(0, 1, 26)
    for Np in (14, 16):
        F = (Np - 2) // 2
        xs, xc = pooled(Np, F, "h", "sampled"), pooled(Np, F, "h", "chord")
        if not len(xs):
            continue
        ax.hist(1 / (1 + xs), bins=bins, density=True, histtype="step", color=cols[Np], lw=1.2, label=f"N'={Np} first order")
        ax.hist(1 / (1 + xc), bins=bins, density=True, histtype="step", color=cols[Np], lw=0.9, ls="--")
        le = exact_decoder(Np, -2)
        if len(le):
            ax.hist(le, bins=bins, density=True, histtype="step", color=cols[Np], lw=0.8, ls=":")
    ax.plot([], [], "k--", lw=0.9, label="chord (first order)")
    ax.plot([], [], "k:", lw=0.8, label="exact uplift")
    ax.set_xlabel(r"decoder eigenvalue $\lambda$")
    ax.set_ylabel("density")
    ax.legend(fontsize=5.8, frameon=False, loc="upper center")
    ax.set_title("(b) first-order decoder", fontsize=8.5, loc="left")
    ax = axs[2]
    for (lp, mk, lab) in ((-2, "o", "resonant"), (0, "s", "central")):
        rs = [r for r in rows if r["lp"] == lp and r["fam"] == "h"]
        if not rs:
            continue
        Nps = [r["Np"] for r in rs]
        ax.plot(Nps, [100 * r["W_total"] for r in rs], mk + "-", color="#d95f02" if lp == -2 else "0.4", ms=4, lw=1,
                label=f"W, {lab}")
        ax.plot(Nps, [100 * r["T2_total"] for r in rs], mk + "--", color="#d95f02" if lp == -2 else "0.4", ms=4, lw=1,
                mfc="none", label=fr"$T_2$, {lab}")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel(r"$N'$")
    ax.set_ylabel("chord / exact first order − 1 (%)")
    ax.legend(fontsize=5.8, frameon=False)
    ax.set_title("(c) decoder statistics", fontsize=8.5, loc="left")
    fig.tight_layout(pad=0.4, w_pad=0.6)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"fig_layer1.{ext}"), dpi=200)
    plt.close(fig)


def main():
    recs = load()
    rows = table(recs)
    print(" N' l' fam seeds  b    e   Emin/E | W: sampled wishart chord (planar, continuum, total %) | T2: sampled wishart chord"
          " (planar, continuum, total %) | P*sqrtK(100) sampled/wishart/chord | Hill(K>100) sampled/chord")
    for r in rows:
        print(f"{r['Np']:3d} {r['lp']:+d} {r['fam']} {r['seeds']} {r['b']:5d} {r['e']:5d} {r['Emin_rel']:.1e} | "
              f"{r['W_sampled']:.4f} {r['W_wishart']:.4f} {r['W_chord']:.4f} ({100*r['W_planar']:+.1f}, {100*r['W_continuum']:+.1f},"
              f" {100*r['W_total']:+.1f}) | {r['T2_sampled']:.4f} {r['T2_wishart']:.4f} {r['T2_chord']:.4f}"
              f" ({100*r['T2_planar']:+.1f}, {100*r['T2_continuum']:+.1f}, {100*r['T2_total']:+.1f}) | "
              f"{r.get('Psqrt100_sampled', float('nan')):.3f}/{r.get('Psqrt100_wishart', float('nan')):.3f}/"
              f"{r.get('Psqrt100_chord', float('nan')):.3f} | {r.get('hill100_sampled', float('nan')):.3f}"
              f" (n={r.get('hill100_n_sampled', 0)})/{r.get('hill100_chord', float('nan')):.3f}")
    with open(os.path.join(DATA, "summary.json"), "w") as fh:
        json.dump(rows, fh, indent=1)
    figure(rows)
    print("wrote", os.path.join(DATA, "summary.json"), "and", os.path.join(FIGS, "fig_layer1.pdf"))


if __name__ == "__main__":
    main()
