#!/usr/bin/env python3
r"""
What localization looks like: decoder angle density vs an isotropic deformation at (nearly) equal T_2, with the
first-order uplift prediction from the old model (progress report 2026-10-13, docs/derivations.md D6 §6).

y = cos^2(2 theta) for the decoder (theta = principal angle between the BPS space and the n_e = 0 subspace);
for the isotropic OU pair (P, P_x) the analogous variable is c^2 = cos^2 of the principal angle between the two BPS
spaces, and T_2m = <y^m> (resp. <c^{2m}>) in both cases.  The OU lag closest to the matched x_* is used (its T_2 is
quoted; it is not exactly equal to the decoder's).  First order: y_1 = ((1 - s^2)/(1 + s^2))^2 over all seeds of
results/data/uplift_2026-10-12/walls_N{N'}.jsonl.

Usage: .venv/bin/python scripts/plot_angle_densities.py --Np 14 15
"""
import argparse
import json
import os
import sys

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from analyze_ou_fourarc import analyze, load  # noqa: E402

WALLS = os.path.join(ROOT, "results/data/uplift_2026-10-12")


def densities(N, nbins_out=25):
    rows = load(N)
    xstar = analyze(N, nboot=20)["decoder"][1]["x_star"]
    lags = rows[0]["ou"]
    L = min(lags, key=lambda o: abs(o["x"] - xstar))["lag"]
    h_ou = np.sum([np.array(o["hist"], float) for r in rows for o in r["ou"] if o["lag"] == L], axis=0)
    m_ou = np.mean([o["moments"][:2] for r in rows for o in r["ou"] if o["lag"] == L], axis=0)
    xL = [o["x"] for o in lags if o["lag"] == L][0]
    h_dec = np.sum([np.array(p["hist"], float) for r in rows for p in r["probes"] if p["s"] == 1], axis=0)
    m_dec = np.mean([p["moments"][:2] for r in rows for p in r["probes"] if p["s"] == 1], axis=0)
    y1 = []
    for line in open(os.path.join(WALLS, f"walls_N{N}.jsonl")):
        lam = np.array(json.loads(line)["first_lam"])
        y1.append((2 * lam - 1) ** 2)
    y1 = np.concatenate(y1)
    nb = len(h_dec)
    f = nb // nbins_out                                       # coarse-grain the stored 200-bin histograms
    edges = np.linspace(0, 1, nbins_out + 1)

    def coarse(h):
        hc = h.reshape(nbins_out, f).sum(1)
        return hc / hc.sum() * nbins_out
    h1, _ = np.histogram(y1, bins=edges)
    out = dict(N=N, edges=edges, dec=coarse(h_dec), ou=coarse(h_ou), first=h1 / h1.sum() * nbins_out,
               T2_dec=m_dec[0], H_dec=m_dec[0] - m_dec[1], T2_ou=m_ou[0], H_ou=m_ou[0] - m_ou[1], x_lag=xL,
               x_star=xstar, T2_first=float(y1.mean()), H_first=float(np.mean(y1 * (1 - y1))))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Np", type=int, nargs="+", default=[14, 15])
    ap.add_argument("--out", default=os.path.join(ROOT, "results/figures/uplift_2026-10-12/fig_angle_densities.pdf"))
    ap.add_argument("--sym", default="H", help="symbol for <y(1-y)> in legends (the 2026-10-13 report uses W)")
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.labelsize": 10, "legend.fontsize": 8.2})
    fig, axes = plt.subplots(1, len(a.Np), figsize=(4.8 * len(a.Np), 3.7), squeeze=False)
    for ax, N in zip(axes[0], a.Np):
        d = densities(N)
        e = d["edges"]
        ax.stairs(d["ou"], e, color="tab:gray", lw=1.4, ls="--",
                  label=rf"isotropic OU, $x={d['x_lag']:.3f}$: $T_2$={d['T2_ou']:.3f}, ${a.sym}$={d['H_ou']:.3f}")
        ax.stairs(d["dec"], e, color="k", lw=1.6,
                  label=rf"decoder: $T_2$={d['T2_dec']:.3f}, ${a.sym}$={d['H_dec']:.3f}")
        ax.stairs(d["first"], e, color="tab:red", lw=1.3, ls=":",
                  label=rf"first-order uplift: ${a.sym}$={d['H_first']:.3f}")
        ax.axvline(0.5, color="k", lw=0.5, alpha=0.4)
        ax.set_yscale("log")
        ax.set_ylim(0.2, 12)
        from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
        ax.yaxis.set_major_locator(FixedLocator([0.3, 0.5, 1, 2, 5, 10]))
        ax.yaxis.set_major_formatter(FixedFormatter(["0.3", "0.5", "1", "2", "5", "10"]))
        ax.yaxis.set_minor_locator(NullLocator())
        ax.set_xlim(0, 1)
        ax.set_xlabel(r"$y=\cos^2 2\theta$  (OU: $\cos^2\theta$)")
        ax.set_ylabel("density")
        ax.set_title(rf"$N'={N}$ (decoder $x_*={d['x_star']:.3f}$)", fontsize=10)
        ax.legend(frameon=False, loc="upper center")
        print(N, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in d.items() if not isinstance(v, np.ndarray)})
    fig.tight_layout()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.savefig(a.out)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
