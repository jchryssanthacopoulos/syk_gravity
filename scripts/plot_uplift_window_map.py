#!/usr/bin/env python3
r"""
Jump map of the uplift anywhere in the fortuity window (docs/derivations.md D10), in the style of CV Figs. 1-3.

fig_jump_map_q{q}: BPS dimensions in the (N, l = 2p - N) plane, CV conventions.
- Light-blue window |l| <= q; exceptional states as white squares; dashed lines bound the restricted window.
- Right margin: the regime of the uplift INTO each row l' (D10 §1).
  - central (|l'| <= q-2): both dressing targets lie outside the window (gapped dressing);
  - resonant (|l'| = q-1): one dressing target is the window edge, the gapless KW multiplet;
  - edge (|l'| = q): the dressing target lies inside the window, delta obstructs, and the radical gives the
    exceptional states.
- A red cross in a cell: the connecting map delta of the uplift into that sector is nonzero.

fig_jump_panels_q3: four uplifts N -> N' = N + 1, one per regime, in columns N and N'.
- Jumps: the n_e = 0 slot comes from l'+1 (down step) and the n_e = 1 slot from l'-1 (up step).
- First-order dressing loops, drawn left of column N:
  - E: n_e = 0 slot -> l' + 2q - 1, then Q0^+ back to the n_e = 1 slot;
  - E^dag: n_e = 1 slot -> l' - 2q + 1, then (Q0^dag)^+ back.
- Loop colour: grey gapped, orange window edge, red inside the window.
- Boxes give rank delta, the admissible dimensions and the edge gap, read from
  results/data/uplift_window_2026-10-26/sectors_p3.jsonl (seed 9000) when present.
Dimensions: q = 3 from uplift.bps_dim_q3 (strip walks + KW exceptional counts, N <= 17).
q = 5: strip walks + exceptional states at N = 7, 11 (KW §6.5, checked to N = 17).

Usage: .venv/bin/python scripts/plot_uplift_window_map.py [--q 3] [--Nmax 16]
"""
import argparse
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from matplotlib.patches import FancyArrowPatch, Rectangle  # noqa: E402

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from uplift import bps_dim_q3, strip_walk_count  # noqa: E402

EXC_Q5 = {7: 1, 11: 1}
WINDOW, OUTSIDE = "#d9e2f7", "#efefef"
REG_COL = {"gapped": "#6b6b6b", "resonant": "#e08214", "obstructed": "#c51b1b"}
plt.rcParams.update({"font.size": 9, "font.family": "serif", "mathtext.fontset": "cm"})


def dim(N, l, q):
    """BPS dimension at (N, l); None where the parity is wrong."""
    if (N + l) % 2 or abs(l) > N:
        return None
    if abs(l) > q:
        return 0
    F = (N + l) // 2
    if q == 3:
        return bps_dim_q3(N, F)
    if abs(l) == q:
        return EXC_Q5.get(N, 0)
    return strip_walk_count(N, q, F)


def obstructed(Np, lp, q):
    """delta != 0 for the uplift into (N', l') (D10 §3).
    - |l'| = q: whenever the slot at -+(q-1) is nonempty.
    - |l'| = q-1 or q+1: when exceptional states of N = N'-1 sit at the far edge, resp. in the slot itself."""
    N = Np - 1
    if abs(lp) == q:
        return (dim(N, lp + 1 if lp < 0 else lp - 1, q) or 0) > 0
    if abs(lp) in (q - 1, q + 1):
        edge = q if lp < 0 else -q
        if abs(lp) == q + 1:
            edge = -q if lp < 0 else q
        return (dim(N, edge, q) or 0) > 0 and (dim(N, lp + 1 if lp < 0 else lp - 1, q) or 0) > 0
    return False


def cell(ax, x, l, d, norm, cmap, q, fs=7.0, w=1.0, edge=None, lw=1.6):
    if d:
        if abs(l) == q:
            ax.add_patch(Rectangle((x - w / 2, l - 0.5), w, 1, fc="white", ec="0.2", lw=0.9, zorder=2))
            ax.text(x, l, f"{d}", ha="center", va="center", fontsize=fs, zorder=3)
        else:
            c = cmap(norm(d))
            ax.add_patch(Rectangle((x - w / 2, l - 0.5), w, 1, fc=c, ec="none", zorder=2))
            lum = 0.299 * c[0] + 0.587 * c[1] + 0.114 * c[2]
            ax.text(x, l, f"{d}", ha="center", va="center", fontsize=fs if d < 1000 else fs - 1.2,
                    color="white" if lum < 0.5 else "black", zorder=3)
    if edge:
        ax.add_patch(Rectangle((x - w / 2, l - 0.5), w, 1, fc="none", ec=edge, lw=lw, zorder=4))


def bands(ax, x0, x1, q, lr):
    ax.add_patch(Rectangle((x0, -q - 0.5), x1 - x0, 2 * q + 1, fc=WINDOW, ec="none", zorder=0))
    ax.add_patch(Rectangle((x0, q + 0.5), x1 - x0, lr - q, fc=OUTSIDE, ec="none", zorder=0))
    ax.add_patch(Rectangle((x0, -lr - 0.5), x1 - x0, lr - q, fc=OUTSIDE, ec="none", zorder=0))
    for y in (q - 0.5, -q + 0.5):
        ax.plot([x0, x1], [y, y], color="0.35", ls=(0, (3, 2)), lw=0.7, zorder=1)


def global_map(fig, ax, q, Nmax, norm, cmap):
    lr = q + 1
    bands(ax, 0.5, Nmax + 0.5, q, lr)
    for N in range(1, Nmax + 1):
        for l in range(-lr, lr + 1):
            d = dim(N, l, q)
            if d is None:
                continue
            cell(ax, N, l, d, norm, cmap, q)
            if N >= 2 and abs(l) <= q + 1 and obstructed(N, l, q):
                ax.plot([N + 0.3], [l + 0.3], marker="x", color=REG_COL["obstructed"], ms=4.2, mew=1.3, zorder=5)
    ax.set_xlim(0.5, Nmax + 0.5)
    ax.set_ylim(-lr - 0.5, lr + 0.5)
    ax.set_xticks(range(1, Nmax + 1))
    ax.set_yticks(range(-lr, lr + 1))
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"$\ell = 2p - N$")
    ax.tick_params(length=2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    xb = Nmax + 0.75
    rows = [((-(q - 2), q - 2), "central", "gapped dressing", REG_COL["gapped"])]
    for sgn in (1, -1):
        rows.append(((sgn * (q - 1), sgn * (q - 1)), "resonant", "via edge multiplet", REG_COL["resonant"]))
        rows.append(((sgn * q, sgn * q), "edge", r"$\delta\neq0$, radical", REG_COL["obstructed"]))
    for (a, b), t1, t2, col in rows:
        lo, hi = min(a, b), max(a, b)
        ax.plot([xb, xb], [lo - 0.42, hi + 0.42], color=col, lw=3, clip_on=False, solid_capstyle="butt")
        ax.text(xb + 0.25, (lo + hi) / 2, f"{t1}: {t2}", color=col, va="center", ha="left", fontsize=7.3, clip_on=False)
    ax.text(xb + 0.25, lr, "outside", color="0.45", va="center", ha="left", fontsize=7.3, clip_on=False)
    ax.text(xb + 0.25, -lr, "outside", color="0.45", va="center", ha="left", fontsize=7.3, clip_on=False)


def load_records(path):
    recs = {}
    if os.path.exists(path):
        for line in open(path):
            r = json.loads(line)
            if r.get("task") == "sectors" and r.get("p") == 3 and r.get("seed") == 9000:
                recs[(r["Np"], r["lp"])] = r
    return recs


def loop(ax, x, y0, y1, col, rad, ls="-", lw=1.4, label=None, lx=None):
    ax.add_patch(FancyArrowPatch((x, y0), (x, y1), connectionstyle=f"arc3,rad={rad}", arrowstyle="-|>",
                                 mutation_scale=9, color=col, lw=lw, ls=ls, zorder=6, shrinkA=3, shrinkB=3))
    if label:
        ax.text(lx, (y0 + y1) / 2, label, color=col, fontsize=8, ha="center", va="center", zorder=7,
                bbox=dict(fc="white", ec="none", pad=0.4, alpha=0.85))


def target_regime(t, q):
    return "gapped" if abs(t) > q else ("resonant" if abs(t) == q else "obstructed")


def panel(ax, N, lp, q, norm, cmap, recs, title, notes=""):
    Np, lr = N + 1, 2 * q + 1
    xN, xP = 0.0, 1.55
    bands(ax, -2.15, 2.25, q, lr)
    for x, n in ((xN, N), (xP, Np)):
        for l in range(-lr, lr + 1):
            d = dim(n, l, q)
            if d is None:
                continue
            cell(ax, x, l, d, norm, cmap, q, fs=7.5, edge="black" if (n == Np and l == lp) else None, lw=2.0)
    for ls_, lab in ((lp + 1, r"$n_e{=}0$"), (lp - 1, r"$n_e{=}1$")):
        if (dim(N, ls_, q) or 0) > 0:
            ax.add_patch(FancyArrowPatch((xN + 0.5, ls_), (xP - 0.5, lp), arrowstyle="-|>", mutation_scale=10,
                                         color="0.1", lw=1.3, zorder=6, shrinkA=1, shrinkB=1))
            ax.text((xN + xP) / 2, (ls_ + lp) / 2 + (0.5 if ls_ > lp else -0.5), lab, fontsize=6.8,
                    ha="center", va="center", zorder=7)
    ta, tb = lp + 1 + 2 * (q - 1), lp - 1 - 2 * (q - 1)
    for src, tgt, xe, xr, rad, lE, lR in ((lp + 1, ta, xN - 0.62, xN - 0.58, -0.5, r"$E$", r"$Q_0^{+}$"),
                                          (lp - 1, tb, xN - 0.62, xN - 0.70, 0.5, r"$E^\dagger$", r"$(Q_0^\dagger)^{+}$")):
        if (dim(N, src, q) or 0) == 0 or abs(tgt) > lr:
            continue
        back = lp - 1 if src == lp + 1 else lp + 1
        reg = target_regime(tgt, q)
        col = REG_COL[reg]
        ax.add_patch(Rectangle((xN - 0.5, tgt - 0.5), 1, 1, fc="none", ec=col, lw=1.2, ls=(0, (2, 1.5)), zorder=5))
        if (dim(N, tgt, q) or 0) == 0:
            ax.text(xN, tgt, {"gapped": "gapped", "resonant": "gapless edge"}.get(reg, ""), fontsize=6.5, color=col,
                    ha="center", va="center", zorder=6)
        else:
            ax.text(-2.1, tgt, {"resonant": "exceptional\n+ gapless edge", "obstructed": "BPS target"}.get(reg, ""),
                    fontsize=6.2, color=col, ha="left", va="center", zorder=6, linespacing=0.95)
        loop(ax, xe, src, tgt, col, rad=rad, label=lE, lx=xN - 1.62)
        ax.add_patch(FancyArrowPatch((xr, tgt), (xr, back), arrowstyle="-|>", mutation_scale=8, color=col,
                                     lw=1.0, ls=(0, (3, 2)), zorder=6, shrinkA=4, shrinkB=4))
        ax.text(xr - 0.08, (tgt + back) / 2 + (0.9 if src == lp + 1 else -0.9), lR, color=col, fontsize=7.3,
                ha="right", va="center", zorder=7, bbox=dict(fc="white", ec="none", pad=0.3, alpha=0.85))
    ax.set_xlim(-2.15, 2.25)
    ax.set_ylim(-lr - 0.5, lr + 0.5)
    ax.set_xticks([xN, xP])
    ax.set_xticklabels([rf"$N={N}$", rf"$N'={Np}$"])
    ax.set_yticks(range(-lr + 1, lr, 2))
    ax.set_ylabel(r"$\ell$", labelpad=1)
    ax.tick_params(length=2, labelsize=7.5)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.set_title(title, fontsize=8.6, loc="left")
    r = recs.get((Np, lp))
    if r:
        rk = r["obstruction"]["h"]["rank"], r["obstruction"]["g"]["rank"]
        txt = rf"rank $\delta$ = {rk[0]} + {rk[1]}" + "\n" + rf"$d = d_h + d_g$ = {r['d_h']} + {r['d_g']}"
        inn = r.get("inner", {}).get("h", {})
        if abs(lp) == q - 1 and inn.get("s2_min"):
            txt += "\n" + rf"edge gap $E_{{\min}}/\bar E$ = {inn['s2_min'] / inn['s2_mean']:.0e}"
        if notes:
            txt += "\n" + notes
        ax.text(2.2, lr + 0.35, txt, fontsize=6.9, va="top", ha="right", zorder=8,
                bbox=dict(fc="white", ec="0.65", lw=0.5, pad=2))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--q", type=int, default=3)
    ap.add_argument("--Nmax", type=int, default=16)
    ap.add_argument("--out", default=os.path.join(ROOT, "results/figures/uplift_window_2026-10-26"))
    ap.add_argument("--data", default=os.path.join(ROOT, "results/data/uplift_window_2026-10-26/sectors_p3.jsonl"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    q = a.q
    cmap = plt.get_cmap("inferno").copy()
    dmax = max((dim(N, l, q) or 0) for N in range(1, a.Nmax + 1) for l in range(-q, q + 1))
    norm = LogNorm(vmin=0.55, vmax=dmax * 2.5)

    fig, ax = plt.subplots(figsize=(7.2, 3.0 if q == 3 else 3.9))
    global_map(fig, ax, q, a.Nmax, norm, cmap)
    fig.subplots_adjust(left=0.07, right=0.765, top=0.97, bottom=0.14)
    fn = os.path.join(a.out, f"fig_jump_map_q{q}")
    for ext in ("pdf", "png"):
        fig.savefig(f"{fn}.{ext}", dpi=220)
    print("wrote", fn)

    if q == 3:
        recs = load_records(a.data)
        fig, axs = plt.subplots(2, 2, figsize=(7.2, 7.4))
        spec = [(11, 0, r"(a) central, $\ell'=0$", ""),
                (11, -2, r"(b) resonant, $\ell'=-2$", ""),
                (9, -2, r"(c) resonant + exchange, $\ell'=-2$", "3 exceptional b-states\nreplace 3 obstructed a-states"),
                (12, -3, r"(d) edge, $\ell'=-3$", r"radical of the antisymmetric" + "\n" + r"form on $H(12,-2)$")]
        for ax, (N, lp, t, nt) in zip(axs.flat, spec):
            panel(ax, N, lp, q, norm, cmap, recs, t, nt)
        handles = [plt.Line2D([], [], color="0.1", lw=1.3, label="jump (walk step)"),
                   plt.Line2D([], [], color=REG_COL["gapped"], lw=1.4, label="dressing, target outside window"),
                   plt.Line2D([], [], color=REG_COL["resonant"], lw=1.4, label="dressing, target = window edge"),
                   plt.Line2D([], [], color=REG_COL["obstructed"], lw=1.4, label="dressing, target inside window"),
                   plt.Line2D([], [], color="0.3", lw=1.1, ls=(0, (3, 2)), label=r"return $Q_0^{+}$ / $(Q_0^\dagger)^{+}$")]
        fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=7.4, frameon=False, bbox_to_anchor=(0.5, 0.0))
        fig.subplots_adjust(left=0.07, right=0.99, top=0.965, bottom=0.085, hspace=0.22, wspace=0.18)
        fn = os.path.join(a.out, "fig_jump_panels_q3")
        for ext in ("pdf", "png"):
            fig.savefig(f"{fn}.{ext}", dpi=220)
        print("wrote", fn)


if __name__ == "__main__":
    main()
