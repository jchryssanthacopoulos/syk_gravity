#!/usr/bin/env python3
r"""
Analyze the chord-invariant campaign (results/data/chord_invariants_2026-09-28/raw_*.jsonl).

Produces (never overwrites raw data):
  results/data/chord_invariants_2026-09-28/summary.csv          per (q, N', model): mean +- s.e. of all invariants
  results/data/chord_invariants_2026-09-28/summary_words.csv    chord-rule test by crossing number
  results/data/chord_invariants_2026-09-28/fits.json            scaling fits
  results/figures/chord_invariants_2026-09-28/*.png, *.pdf      figures
and prints markdown tables used in research/notes/chord_invariants_2026-09-28.md.

Usage: .venv/bin/python scripts/analyze_chord_invariants.py [--data DIR] [--figs DIR]
"""
import argparse
import glob
import json
import os
import sys
from collections import defaultdict
from math import comb, pi, sqrt

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from chord_invariants import free_compression_moments  # noqa: E402

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"     # validated categorical slots 1-3 (dataviz palette)
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


# ----------------------------------------------------------------------------- loading / reduction
def load(data):
    rows = []
    for f in sorted(glob.glob(os.path.join(data, "raw_q*_N*.jsonl"))):
        with open(f) as fh:
            rows += [json.loads(line) for line in fh if line.strip()]
    return rows


def per_realization(row):
    """Reduce one realization to scalars (averages over modes / pairs / triples)."""
    out = dict(q=row["q"], Np=row["Np"], model=row["model"], a=row["a"], b=row["b"], d=row["d"], D=row["D"])
    S = row["single"]
    rho = np.array([s["rho"] for s in S])                    # modes x 6 (index n-1)
    out["r"] = float(np.mean(rho[:, 1]))
    for n in range(2, 7):
        vals = rho[:, n - 1]
        out[f"rho{n}"] = float(np.nanmean(vals))
    dec = [s for s in S if s["mode"] == row["Np"] - 1][0]
    out["r_dec"] = dec["r"]
    out["atoms_frac"] = float(np.mean([(s["n0"] + s["n1"]) / row["d"] for s in S]))
    if row["pairs"]:
        for key in ("c", "F", "X"):
            out[key] = float(np.mean([p[key] for p in row["pairs"]]))
        out["XoverF"] = float(np.mean([p["X"] / p["F"] for p in row["pairs"]]))
    # chord rule: average ratio per crossing number, normalized by X^cr (using this realization's X)
    by_cr = defaultdict(list)
    for t in row["triples"]:
        for w in t["words"]:
            by_cr[w["cr"]].append(w["ratio"])
    for cr, v in by_cr.items():
        out[f"w{cr}"] = float(np.mean(v))
        out[f"w{cr}_spread"] = float(np.std(v))
    for k, mm in row.get("multimode", {}).items():
        out[f"r_k{k}"] = mm["single"]["r"]
        out[f"rho4_k{k}"] = mm["single"]["rho"][3]
        out[f"X_k{k}"] = mm["pair"]["X"]
        out[f"F_k{k}"] = mm["pair"]["F"]
        out[f"b_k{k}"] = mm["single"]["b_uv"]
    return out


def aggregate(reals):
    groups = defaultdict(list)
    for r in reals:
        groups[(r["q"], r["Np"], r["model"])].append(r)
    table = {}
    for key, rs in sorted(groups.items()):
        keys = sorted(set().union(*[set(r) for r in rs]) - {"q", "Np", "model"})
        agg = dict(q=key[0], Np=key[1], model=key[2], n=len(rs))
        for k in keys:
            v = np.array([r[k] for r in rs if k in r and r[k] is not None], float)
            v = v[np.isfinite(v)]
            if len(v) == 0:
                continue
            agg[k] = float(v.mean())
            agg[k + "_se"] = float(v.std(ddof=1) / np.sqrt(len(v))) if len(v) > 1 else float("nan")
        table[key] = agg
    return table


# ----------------------------------------------------------------------------- Wachter law with parameter t
def wachter_cdf(t, b, xs):
    """CDF (per-BPS normalization, including atoms) of the free compression of Bern(b) by a projection of trace t."""
    lm = (sqrt(t * (1 - b)) - sqrt(b * (1 - t))) ** 2
    lp = (sqrt(t * (1 - b)) + sqrt(b * (1 - t))) ** 2
    p0, p1 = max(0.0, 1 - b / t), max(0.0, 1 - (1 - b) / t)
    g = np.linspace(lm, lp, 20001)[1:-1]
    f = np.sqrt(np.clip((lp - g) * (g - lm), 0, None)) / (2 * pi * t * g * (1 - g))
    c = np.concatenate([[0.0], np.cumsum((f[1:] + f[:-1]) / 2 * np.diff(g))])
    c = c / c[-1] * (1 - p0 - p1)                             # continuous mass
    out = np.interp(xs, g, c, left=0.0, right=1 - p0 - p1) + p0 * (np.asarray(xs) >= 0)
    return np.where(np.asarray(xs) >= 1, 1.0, out)


def wachter_density(t, b, xs):
    lm = (sqrt(t * (1 - b)) - sqrt(b * (1 - t))) ** 2
    lp = (sqrt(t * (1 - b)) + sqrt(b * (1 - t))) ** 2
    xs = np.asarray(xs)
    f = np.zeros_like(xs)
    m = (xs > lm) & (xs < lp)
    f[m] = np.sqrt((lp - xs[m]) * (xs[m] - lm)) / (2 * pi * t * xs[m] * (1 - xs[m]))
    return f


def ks(lam, t, b, tol=1e-9):
    """KS distance between the full empirical law (atoms included) and Wachter(t, b)."""
    lam = np.sort(np.where(lam > 1 - tol, 1.0, np.where(lam < tol, 0.0, lam)))
    n = len(lam)
    F = wachter_cdf(t, b, lam)
    # empirical CDF just after / just before each point; at atoms compare both sides of the jump
    emp_hi = np.searchsorted(lam, lam, side="right") / n
    emp_lo = np.searchsorted(lam, lam, side="left") / n
    F_lo = wachter_cdf(t, b, lam - 1e-12)
    return float(max(np.max(np.abs(emp_hi - F)), np.max(np.abs(emp_lo - F_lo))))


def ks_interior(lam, t, b, tol=1e-9):
    """KS distance of the interior (0<lam<1) empirical law vs the continuous part of Wachter(t, b), both unit mass."""
    x = np.sort(lam[(lam > tol) & (lam < 1 - tol)])
    p0, p1 = max(0.0, 1 - b / t), max(0.0, 1 - (1 - b) / t)
    cont = 1 - p0 - p1
    F = (wachter_cdf(t, b, x) - p0) / cont
    n = len(x)
    return float(max(np.max(np.abs(np.arange(1, n + 1) / n - F)), np.max(np.abs(np.arange(0, n) / n - F))))


def uv_refs(Np, P, T1, T2, T3=None):
    """Commuting (unprojected) reference values on the full charge-P sector: X_UV, F_UV, and (if T3) w3_UV."""
    from q_scan import subset_index
    masks, _ = subset_index(Np, P)
    def slot(T):
        x = np.ones(len(masks))
        for i in T:
            x *= np.array([1.0 - ((m >> i) & 1) for m in masks])
        return x - x.mean()
    y1, y2 = slot(T1), slot(T2)
    v1, v2 = np.mean(y1 ** 2), np.mean(y2 ** 2)
    out = dict(X=float(np.mean(y1 * y2 * y1 * y2) / (v1 * v2)), F=float(np.mean(y1 * y1 * y2 * y2) / (v1 * v2)))
    if T3 is not None:
        y3 = slot(T3)
        out["w3"] = float(np.mean(y1 * y2 * y3 * y1 * y2 * y3) / (v1 * v2 * np.mean(y3 ** 2)))
    return out


# ----------------------------------------------------------------------------- fits
def fit_line(x, y, w=None):
    x, y = np.asarray(x, float), np.asarray(y, float)
    A = np.vstack([np.ones_like(x), x]).T
    W = np.ones_like(x) if w is None else np.asarray(w, float)
    coef, *_ = np.linalg.lstsq(A * W[:, None], y * W, rcond=None)
    res = y - A @ coef
    return coef, float(np.sum((res * W) ** 2))


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=os.path.join(ROOT, "results", "data", "chord_invariants_2026-09-28"))
    ap.add_argument("--figs", default=os.path.join(ROOT, "results", "figures", "chord_invariants_2026-09-28"))
    args = ap.parse_args()
    os.makedirs(args.figs, exist_ok=True)
    rows = load(args.data)
    reals = [per_realization(r) for r in rows]
    T = aggregate(reals)

    # ---------------- summary.csv
    cols = ["q", "Np", "model", "n", "d", "a", "b", "r", "r_se", "rho4", "rho4_se", "rho6", "rho6_se",
            "X", "X_se", "F", "F_se", "c", "c_se", "atoms_frac", "w0", "w1", "w2", "w3",
            "r_k2", "X_k2", "r_k3", "X_k3"]
    with open(os.path.join(args.data, "summary.csv"), "w") as f:
        f.write(",".join(cols) + "\n")
        for key, g in T.items():
            f.write(",".join(str(g.get(c, "")) for c in cols) + "\n")

    def fmt(g, k, p=3):
        v, e = g.get(k), g.get(k + "_se")
        if v is None:
            return "–"
        return f"{v:.{p}f}" + (f"±{e:.{p}f}" if e is not None and np.isfinite(e) else "")

    print("\n### Table 1: single-operator and pair invariants (mean ± s.e. over realizations)\n")
    print("| q | N′ | model | n | a | r | ρ₄ (a³) | ρ₆ (a⁵) | X | F | c |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for (q, Np, model), g in T.items():
        a = g["a"]
        print(f"| {q} | {Np} | {model} | {g['n']} | {a:.3f} | {fmt(g,'r')} | {fmt(g,'rho4')} ({a**3:.3f}) | "
              f"{fmt(g,'rho6')} ({a**5:.3f}) | {fmt(g,'X')} | {fmt(g,'F')} | {fmt(g,'c')} |")

    # ---------------- effective-Wachter test: a_n = rho_n^{1/(n-1)}
    print("\n### Table 2: effective parameters a_n = ρ_n^{1/(n−1)} (Wachter(t) ⇒ a_2 = a_4 = a_6 = t; Haar ⇒ t = a)\n")
    print("| q | N′ | model | a | a₂ = r | a₄ | a₆ | X |")
    print("|---|---|---|---|---|---|---|---|")
    for (q, Np, model), g in T.items():
        a4 = g["rho4"] ** (1 / 3) if g.get("rho4", -1) > 0 else float("nan")
        a6 = g["rho6"] ** (1 / 5) if g.get("rho6", -1) > 0 else float("nan")
        print(f"| {q} | {Np} | {model} | {g['a']:.3f} | {g['r']:.3f} | {a4:.3f} | {a6:.3f} | {g.get('X', float('nan')):.3f} |")

    # ---------------- chord rule (words)
    print("\n### Table 3: chord (q-Gaussian) rule for 6-letter pair words: τ(w)/Πτ(Â²) vs X^cr\n")
    print("cr=3 discriminates mechanisms: free compression (Haar) ⇒ X²;  chord/q-Gaussian rule ⇒ X³.  e₃ ≡ ln w₃ / ln X.\n")
    print("| q | N′ | model | X | cr=0 (pred 1) | cr=1 (pred X) | cr=2 (pred X²) | cr=3 (X³) | X² | e₃ |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    with open(os.path.join(args.data, "summary_words.csv"), "w") as f:
        f.write("q,Np,model,X,w0,w1,w2,w3,w0_se,w1_se,w2_se,w3_se\n")
        for (q, Np, model), g in T.items():
            if "w0" not in g:
                continue
            X = g["X"]
            cells = [f"{g[f'w{c}']:.3f} ({X**c:.3f})" if f"w{c}" in g else "–" for c in range(4)]
            e3 = np.log(g["w3"]) / np.log(X) if "w3" in g and 0 < X < 1 and g["w3"] > 0 else float("nan")
            g["e3"] = e3
            print(f"| {q} | {Np} | {model} | {X:.3f} | " + " | ".join(cells) + f" | {X**2:.3f} | {e3:.2f} |")
            f.write(",".join(str(v) for v in [q, Np, model, X] + [g.get(f"w{c}", "") for c in range(4)]
                             + [g.get(f"w{c}_se", "") for c in range(4)]) + "\n")

    # ---------------- multi-mode probes
    print("\n### Table 4: multi-mode slot probes Π_T, |T| = k (heavier UV operators)\n")
    print("| q | N′ | model | a | r (k=1) | r (k=2) | r (k=3) | X (k=1) | X (k=2) | X (k=3) |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for (q, Np, model), g in T.items():
        print(f"| {q} | {Np} | {model} | {g['a']:.3f} | {fmt(g,'r')} | {fmt(g,'r_k2')} | {fmt(g,'r_k3')} | "
              f"{fmt(g,'X')} | {fmt(g,'X_k2')} | {fmt(g,'X_k3')} |")

    # ---------------- scaling fits (q=3 SYK): y ~ a^s  and  y vs N' (exp vs power)
    fits = {}
    for q in sorted(set(k[0] for k in T)):
        pts = [(Np, g) for (qq, Np, model), g in T.items() if qq == q and model == "syk" and Np >= 7]
        if len(pts) < 3:
            continue
        Nn = np.array([p[0] for p in pts]); A = np.array([p[1]["a"] for p in pts])
        for name in ["r", "X", "rho4", "atoms_frac"]:
            Y = np.array([p[1].get(name, np.nan) for p in pts])
            ok = np.isfinite(Y) & (Y > 0)
            if ok.sum() < 3:
                continue
            (c0, s), rss = fit_line(np.log(A[ok]), np.log(Y[ok]))
            e = dict(power_of_a=float(s), prefactor=float(np.exp(c0)), rss=rss, Np=Nn[ok].tolist())
            one_minus = 1 - Y[ok]
            if np.all(one_minus > 0):
                (c1, s1), rss1 = fit_line(np.log(1 - A[ok]), np.log(one_minus))
                e["one_minus_power_of_one_minus_a"] = float(s1)
            (ce, be), rsse = fit_line(Nn[ok], np.log(Y[ok]))
            (cp, bp), rssp = fit_line(np.log(Nn[ok]), np.log(Y[ok]))
            e.update(exp_rate=float(-be), exp_rss=rsse, power_exponent=float(-bp), power_rss=rssp)
            fits[f"q{q}_{name}"] = e
    with open(os.path.join(args.data, "fits.json"), "w") as f:
        json.dump(fits, f, indent=1)
    print("\n### Table 5: scaling fits (SYK, N′ ≥ 7): y ≈ C·a^s ; also exp-in-N′ vs power-in-N′ RSS\n")
    print("| series | s (y∝a^s; Haar: 1) | C | RSS | exp rate / RSS | power exp. / RSS |")
    print("|---|---|---|---|---|---|")
    for k, e in fits.items():
        print(f"| {k} | {e['power_of_a']:.3f} | {e['prefactor']:.3f} | {e['rss']:.1e} | "
              f"{e['exp_rate']:.3f} / {e['exp_rss']:.1e} | {e['power_exponent']:.3f} / {e['power_rss']:.1e} |")

    # ---------------- spectra: KS vs Wachter(a) and Wachter(r), and moments vs Wachter(r)
    print("\n### Table 6: decoder spectrum (mode N′−1) vs Wachter(a) and Wachter(t=r), pooled over realizations\n")
    print("| q | N′ | model | a | r | KS W(a) | KS W(r) | KS_int W(a) | KS_int W(r) | max rel. mom. err W(a) | W(r) |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    spec_dir = os.path.join(args.data, "spectra")
    spec_store = {}
    for (q, Np, model), g in T.items():
        files = sorted(glob.glob(os.path.join(spec_dir, f"q{q}_N{Np}_{model}_s*.npy")))
        if not files:
            continue
        lam = np.concatenate([np.load(f) for f in files])
        b = 1 - (Np // 2) / Np
        t = float(np.clip(g["r_dec"], 1e-6, 1 - 1e-6))
        mw = free_compression_moments(b, t, 6)
        ma = free_compression_moments(b, g["a"], 6)
        mm = [float(np.mean(lam ** k)) for k in range(1, 7)]
        err = max(abs(mm[k] - mw[k]) / mw[k] for k in range(6))
        erra = max(abs(mm[k] - ma[k]) / ma[k] for k in range(6))
        spec_store[(q, Np, model)] = (lam, g["a"], t, b)
        g.update(ks_a=ks(lam, g["a"], b), ks_r=ks(lam, t, b), ksi_a=ks_interior(lam, g["a"], b),
                 ksi_r=ks_interior(lam, t, b), momerr_a=erra, momerr_r=err)
        print(f"| {q} | {Np} | {model} | {g['a']:.3f} | {t:.3f} | {g['ks_a']:.3f} | {g['ks_r']:.3f} | "
              f"{g['ksi_a']:.3f} | {g['ksi_r']:.3f} | {erra:.3f} | {err:.3f} |")

    # ---------------- rigidity fractions eta = (SYK - Haar)/(UV - Haar): 0 = free/Haar, 1 = commuting UV
    print("\n### Table 7: rigidity fractions η = (SYK − Haar)/(UV − Haar)  (0 = free, 1 = commuting)\n")
    print("| q | N′ | a | η_r (k=1) | η_r (k=2) | η_r (k=3) | η_X (k=1) | η_X (k=2) | η_X (k=3) | X_UV k=1,2,3 |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for (q, Np, model), g in list(T.items()):
        if model != "syk" or (q, Np, "haar") not in T:
            continue
        h = T[(q, Np, "haar")]
        P = Np // 2
        cells_r, cells_X, uvs = [], [], []
        for k in (1, 2, 3):
            rk, Xk = ("r", "X") if k == 1 else (f"r_k{k}", f"X_k{k}")
            if rk not in g or rk not in h or 2 * k > Np or P > Np - k:
                cells_r.append("–"); cells_X.append("–"); continue
            uv = uv_refs(Np, P, list(range(Np - k, Np)) if k > 1 else [Np - 1], list(range(0, k)))
            er = (g[rk] - h[rk]) / (1 - h[rk])
            ex = (g[Xk] - h[Xk]) / (uv["X"] - h[Xk]) if abs(uv["X"] - h[Xk]) > 1e-9 else float("nan")
            g[f"eta_r_k{k}"], g[f"eta_X_k{k}"] = er, ex
            cells_r.append(f"{er:.3f}"); cells_X.append(f"{ex:.3f}"); uvs.append(f"{uv['X']:.3f}")
        print(f"| {q} | {Np} | {g['a']:.3f} | " + " | ".join(cells_r) + " | " + " | ".join(cells_X) + " | " + ", ".join(uvs) + " |")

    make_figures(T, spec_store, args.figs)


# ----------------------------------------------------------------------------- figures
def make_figures(T, spec_store, figs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "legend.frameon": False,
                         "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "savefig.facecolor": "#fcfcfb"})

    def series(q, model, key):
        pts = sorted((Np, g) for (qq, Np, mm), g in T.items() if qq == q and mm == model and key in g)
        a = np.array([g["a"] for _, g in pts]); y = np.array([g[key] for _, g in pts])
        e = np.array([g.get(key + "_se", 0.0) for _, g in pts]); e = np.nan_to_num(e)
        return a, y, e, [Np for Np, _ in pts]

    # Fig 1: r and X vs a
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
    for ax, key, title in [(axs[0], "r", "surviving variance  r = κ₂ / κ₂(UV)"),
                           (axs[1], "X", "crossing weight  X = τ(ÂᵢÂⱼÂᵢÂⱼ)/τ(Âᵢ²)τ(Âⱼ²)")]:
        g = np.linspace(0.25, 1, 50)
        ax.plot(g, g, color=INK2, lw=1.2, ls="--")
        ax.text(0.30, 0.25, "free / Haar: y = a", color=INK2, fontsize=8)
        ax.axhline(1.0, color=INK2, lw=1.0, ls=":")
        ax.text(0.27, 1.01, "commuting (UV): y = 1", color=INK2, fontsize=8, va="bottom")
        for q, col, mk in [(3, BLUE, "o"), (5, ORANGE, "s")]:
            a, y, e, N = series(q, "syk", key)
            if len(a):
                ax.errorbar(a, y, yerr=e, color=col, marker=mk, ms=6, lw=2, capsize=2, label=f"SYK q={q}")
                ax.annotate(f"q={q}", (a[-1], y[-1]), textcoords="offset points", xytext=(6, -12), color=INK, fontsize=8)
        a3, y3, e3, _ = series(3, "haar", key)
        a5, y5, e5, _ = series(5, "haar", key)
        ah, yh = np.concatenate([a3, a5]), np.concatenate([y3, y5])
        ax.plot(ah, yh, ls="none", marker="^", ms=6, color=AQUA, mec="#fcfcfb", mew=1, label="Haar P_B (null model)")
        ax.set_xlabel("BPS fraction  a = d/D"); ax.set_title(title, fontsize=10, color=INK)
        ax.set_xlim(0.25, 1.02)
    axs[0].set_ylabel("invariant")
    axs[1].legend(loc="lower right", fontsize=8)
    fig.suptitle("Decoder chord invariants vs BPS fraction: SYK departs from the free (Haar) line", fontsize=11)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(figs, f"fig1_r_X_vs_a.{ext}"), dpi=160)
    plt.close(fig)

    # Fig 2: effective Wachter parameters a_n vs N' (q=3 SYK)
    pts = sorted((Np, g) for (q, Np, m), g in T.items() if q == 3 and m == "syk")
    if pts:
        fig, ax = plt.subplots(figsize=(6.4, 4.2))
        N = np.array([p[0] for p in pts]); a = np.array([p[1]["a"] for p in pts])
        ax.plot(N, a, color=INK2, lw=1.2, ls="--", label="a = d/D (Wachter/Haar)")
        for key, n, col, mk in [("r", 2, BLUE, "o"), ("rho4", 4, ORANGE, "s"), ("rho6", 6, AQUA, "D")]:
            y = np.array([p[1][key] for p in pts]) if key != "r" else np.array([p[1]["r"] for p in pts])
            an = np.where(y > 0, y, np.nan) ** (1 / (n - 1))
            ax.plot(N, an, color=col, marker=mk, ms=6, lw=2, label=f"a_{n} = ρ_{n}^(1/{n-1})")
        Xs = np.array([p[1].get("X", np.nan) for p in pts])
        ax.plot(N, Xs, color=INK, marker="x", ms=6, lw=1, ls="-.", label="X (crossing weight)")
        ax.set_xlabel("N′ (enlarged size), q = 3"); ax.set_ylabel("effective parameter")
        ax.set_title("Single-mode law ≈ Wachter(t = r), but crossings are not free (X > r)", fontsize=10)
        ax.legend(fontsize=8)
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fig.savefig(os.path.join(figs, f"fig2_effective_parameters_q3.{ext}"), dpi=160)
        plt.close(fig)

    # Fig 3: chord rule test.  (a) w_cr / X^cr for cr = 1, 2 (both mechanisms predict 1);
    #                          (b) e3 = ln w3 / ln X: free compression -> 2, chord/q-Gaussian -> 3
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.0))
    ax = axs[0]
    for cr, col, mk in [(1, BLUE, "o"), (2, ORANGE, "s")]:
        pts = sorted((Np, g) for (q, Np, m), g in T.items() if q == 3 and m == "syk" and f"w{cr}" in g)
        N = [p[0] for p in pts]; y = [p[1][f"w{cr}"] / p[1]["X"] ** cr for p in pts]
        ax.plot(N, y, color=col, marker=mk, ms=6, lw=2, label=f"cr = {cr}:  w/X^{cr}")
    pts = sorted((Np, g) for (q, Np, m), g in T.items() if q == 3 and m == "syk" and "w0" in g)
    ax.plot([p[0] for p in pts], [p[1]["w0"] for p in pts], color=AQUA, marker="D", ms=6, lw=2, label="cr = 0:  w")
    ax.axhline(1, color=INK2, ls="--", lw=1.2)
    ax.set_ylim(0.85, 1.1); ax.set_xlabel("N′ (q = 3, SYK)"); ax.set_ylabel("normalized word value")
    ax.set_title("(a) one or two crossings: multiplicative (1 = exact)", fontsize=10); ax.legend(fontsize=8)
    ax = axs[1]
    for model, col, mk, lab in [("syk", BLUE, "o", "SYK q=3"), ("haar", AQUA, "^", "Haar q=3")]:
        pts = sorted((Np, g) for (q, Np, m), g in T.items() if q == 3 and m == model and np.isfinite(g.get("e3", np.nan)) and Np >= 7)
        ax.plot([p[0] for p in pts], [p[1]["e3"] for p in pts], color=col, marker=mk, ms=6, lw=2, label=lab)
    ax.axhline(2, color=INK2, ls="--", lw=1.2); ax.text(7, 2.03, "free compression: X²", color=INK2, fontsize=8)
    ax.axhline(3, color=INK2, ls=":", lw=1.2); ax.text(7, 3.03, "chord (q-Gaussian) rule: X³", color=INK2, fontsize=8)
    ax.set_ylim(1.6, 3.3); ax.set_xlabel("N′"); ax.set_ylabel("e₃ = ln w₃ / ln X  (word ijkijk)")
    ax.set_title("(b) three mutually crossing chords", fontsize=10); ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(figs, f"fig3_chord_rule_q3.{ext}"), dpi=160)
    plt.close(fig)

    # Fig 4: decoder spectra vs Wachter(a), Wachter(r): full CDFs (atoms visible as jumps)
    sel = [k for k in [(3, 10, "syk"), (3, 12, "syk"), (3, 14, "syk"), (3, 15, "syk")] if k in spec_store]
    if sel:
        fig, axs = plt.subplots(1, len(sel), figsize=(3.4 * len(sel), 3.5), sharey=True)
        axs = np.atleast_1d(axs)
        x = np.linspace(0, 1, 4001)
        for ax, key in zip(axs, sel):
            lam, a, t, b = spec_store[key]
            ls_ = np.sort(lam)
            ax.step(ls_, np.arange(1, len(ls_) + 1) / len(ls_), where="post", color=INK, lw=2, label="SYK (pooled)")
            ax.plot(x, wachter_cdf(a, b, x), color=BLUE, lw=2, ls="--", label="Wachter(a)")
            ax.plot(x, wachter_cdf(t, b, x), color=ORANGE, lw=2, ls="-.", label="Wachter(t = r)")
            ax.set_title(f"N′={key[1]}: a={a:.2f}, r={t:.2f}", fontsize=9)
            ax.set_xlabel("λ (decoder eigenvalue)")
        axs[0].set_ylabel("cumulative fraction of BPS states")
        axs[-1].legend(fontsize=8, loc="lower right")
        fig.suptitle("Decoder spectra (q=3): Wachter(t=r) matches the moments but puts wall mass in atoms; "
                     "SYK spreads it continuously", fontsize=10)
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fig.savefig(os.path.join(figs, f"fig4_spectra_vs_wachter.{ext}"), dpi=160)
        plt.close(fig)

    # Fig 5: rigidity fractions for k-mode probes (q=3)
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.0), sharey=True)
    for ax, key, lab in [(axs[0], "eta_r", "η_r = (r − r_Haar)/(1 − r_Haar)"),
                         (axs[1], "eta_X", "η_X = (X − X_Haar)/(X_UV − X_Haar)")]:
        for k, col, mk in [(1, BLUE, "o"), (2, ORANGE, "s"), (3, AQUA, "D")]:
            kk = f"{key}_k{k}"
            pts = sorted((Np, g) for (q, Np, m), g in T.items() if q == 3 and m == "syk" and kk in g and Np >= 9)
            if not pts:
                continue
            N = [p[0] for p in pts]; y = [p[1][kk] for p in pts]
            ax.plot(N, y, color=col, marker=mk, ms=6, lw=2, label=f"k = {k} modes")
            ax.annotate(f"k={k}", (N[-1], y[-1]), textcoords="offset points", xytext=(6, 0), fontsize=8, color=INK)
        ax.axhline(0, color=INK2, ls="--", lw=1.2); ax.axhline(1, color=INK2, ls=":", lw=1.2)
        ax.set_xlabel("N′ (q = 3, SYK)"); ax.set_title(lab, fontsize=10)
    axs[0].set_ylabel("rigidity (0 = free/Haar, 1 = commuting UV)"); axs[0].legend(fontsize=8)
    fig.suptitle("Heavier slot probes (k modes): surviving variance moves toward free; crossings do not", fontsize=10)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(figs, f"fig5_multimode_rigidity_q3.{ext}"), dpi=160)
    plt.close(fig)

    # Fig 6: matched-a comparison q=3 vs q=5: rigidity eta_r and eta_X vs a
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.0), sharey=True)
    for ax, key in [(axs[0], "eta_r_k1"), (axs[1], "eta_X_k1")]:
        for q, col, mk in [(3, BLUE, "o"), (5, ORANGE, "s")]:
            pts = sorted((g["a"], Np, g) for (qq, Np, m), g in T.items() if qq == q and m == "syk" and key in g)
            if not pts:
                continue
            ax.plot([p[0] for p in pts], [p[2][key] for p in pts], color=col, marker=mk, ms=6, lw=2, label=f"SYK q={q}")
            for A_, N_, g_ in pts:
                ax.annotate(f"{N_}", (A_, g_[key]), textcoords="offset points", xytext=(4, 4), fontsize=7, color=INK2)
        ax.axhline(0, color=INK2, ls="--", lw=1.2)
        ax.set_xlabel("BPS fraction a (labels: N′)"); ax.invert_xaxis()
        ax.set_title("η_r (surviving variance)" if key.startswith("eta_r") else "η_X (crossing weight)", fontsize=10)
    axs[0].set_ylabel("rigidity (0 = free/Haar)"); axs[0].legend(fontsize=8)
    fig.suptitle("Matched-a comparison: single-mode probe, q = 3 vs q = 5", fontsize=10)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(figs, f"fig6_matched_a_q3_q5.{ext}"), dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
