#!/usr/bin/env python3
r"""
Analysis of results/data/uplift_window_2026-10-26 (scripts/run_uplift_window.py; docs/derivations.md D10).

Checks (printed and written to summary.json):
  dims      d_h + d_g = direct BPS dimension = known count (strip walks + exceptional states)
  ranks     connecting-map ranks against the D10 §3 predictions (q = 3 and q = 5):
              |l'| = q-1 or q+1: rank = e(N) (exceptional states of N at the far edge, resp. in the slot);
              |l'| = q: rank = dim H(N, -+(q-1)) - e(N');  otherwise 0
  uplift    residual outside the direct BPS space, rank, closed-form vs direct decoder
  mirror    spec lambda(F) = 1 - spec lambda(N' - F) for the same realization
Tables:
  decoder   per (N', l'): T2 first / self / exact / mean field, H first / exact, T4, T4 - (T2 - H_first), walls
  resonance per (N', l', family): regime, inner E_min/E_mean, first-order and exact dressing-cost tails
Figures (results/figures/uplift_window_2026-10-26/):
  fig_resonance.pdf     (a) inner gap vs N by regime; (b) CCDF of first-order dressing costs, resonant vs gapped,
                        with the chord tail K^{-1/2}; (c) exact vs first-order decoder in a resonant sector
  fig_decoder_sectors.pdf   exact and first-order decoder densities per sector at N' = 13, 14
  fig_T4_sectors.pdf    T4 against T2 - H_first in every sector

Usage: .venv/bin/python scripts/analyze_uplift_window.py
"""
import json
import os
import sys
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from uplift import EXCEPTIONAL_Q3, strip_walk_count  # noqa: E402

DATA = os.path.join(ROOT, "results/data/uplift_window_2026-10-26")
FIGS = os.path.join(ROOT, "results/figures/uplift_window_2026-10-26")
EXC = {3: EXCEPTIONAL_Q3, 5: {7: 1, 11: 1}}
plt.rcParams.update({"font.size": 9, "font.family": "serif", "mathtext.fontset": "cm"})
REG_COL = {"gapped": "#6b6b6b", "resonant": "#e08214", "obstructed": "#c51b1b"}


def load(p):
    fn = os.path.join(DATA, f"sectors_p{p}.jsonl")
    if not os.path.exists(fn):
        return []
    return [json.loads(x) for x in open(fn)]


def chord_tail(Nps=(12, 14, 16), nphi=4000):
    """Parameter-free chord tail of the resonant h family (s_out = -1/6, s_in = 1/2, D10 §4): each inner mode with a large
    jump g contributes one large eigenvalue, so P(sigma^2 > K) ~ r w(g > K) (r = e/b).  Cached in chord_tail.json."""
    fn = os.path.join(DATA, "chord_tail.json")
    cache = json.load(open(fn)) if os.path.exists(fn) else {}
    from chord_dressing import chord_jump_law
    Ks = np.logspace(0, 6, 61)
    for Np in Nps:
        if str(Np) in cache:
            continue
        g, w, P = chord_jump_law(Np, "h", nphi=nphi, F=(Np - 2) // 2)
        cache[str(Np)] = dict(K=Ks.tolist(), P=[float(P["r"] * np.sum(w[g > K])) for K in Ks], r=P["r"], b=P["b"],
                              e=P["e"], lam=P["lam"], s_out=P["s_out"], s_in=P["s_in"])
        json.dump(cache, open(fn, "w"), indent=1)
    return cache


def hill(x, kmin):
    """Hill (maximum-likelihood) estimate of alpha in P(X > K) ~ K^-alpha from the samples above kmin; returns (alpha, n)."""
    x = np.asarray(x, float)
    t = x[x > kmin]
    return (len(t) / np.sum(np.log(t / kmin)), len(t)) if len(t) > 2 else (float("nan"), len(t))


def predicted_ranks(r):
    p, N, lp = r["p"], r["N"], r["lp"]
    e = EXC[p].get(N, 0)
    ep = EXC[p].get(N + 1, 0)
    if lp in (-(p - 1), -(p + 1)):
        return (e, 0)
    if lp in (p - 1, p + 1):
        return (0, e)
    if lp == -p:
        F = (N + 1 - p) // 2                        # old degree of the slot at level 1 - p
        return (strip_walk_count(N, p, F) - ep, 0)
    if lp == p:
        F = (N + p - 1) // 2
        return (0, strip_walk_count(N, p, F) - ep)
    return (0, 0)


def checks(recs):
    out = defaultdict(list)
    for r in recs:
        tag = (r["p"], r["Np"], r["lp"], r["seed"])
        d = r["d_h"] + r["d_g"]
        if r.get("expected_d") is not None and d != r["expected_d"]:
            out["dims_fail"].append(tag)
        if "d_direct" in r and r["d_direct"] != d:
            out["direct_fail"].append(tag)
        rk = (r["obstruction"]["h"]["rank"], r["obstruction"]["g"]["rank"])
        if r["N"] <= 17 and rk != predicted_ranks(r):
            out["rank_fail"].append((tag, rk, predicted_ranks(r)))
        if "uplift_residual" in r:
            out["max_residual"].append(r["uplift_residual"])
            out["max_decoder_dev"].append(r.get("decoder_dev", 0.0))
            if r["uplift_rank"] != d:
                out["uplift_rank_fail"].append(tag)
        if "split_overlap" in r:
            out["max_split_overlap"].append(r["split_overlap"])
        out["n"].append(1)
    summ = {k: (len(v) if k == "n" else (max(v) if k.startswith("max") and v else v)) for k, v in out.items()}
    by = {(r["p"], r["Np"], r["seed"], r["F"]): r for r in recs if r["d_h"] + r["d_g"]}
    mirror = []
    for (p, Np, seed, F), r in by.items():
        m = by.get((p, Np, seed, Np - F))
        if m and F < Np - F:
            a, b = np.array(r["lam"]), np.sort(1 - np.array(m["lam"]))
            mirror.append(float(np.abs(a - b).max()) if len(a) == len(b) else np.inf)
    summ["mirror_pairs"], summ["mirror_max_dev"] = len(mirror), (max(mirror) if mirror else None)
    return summ


def agg(recs, key):
    g = defaultdict(list)
    for r in recs:
        if r["d_h"] + r["d_g"]:
            g[(r["Np"], r["lp"])].append(r)
    return dict(sorted(g.items()))


def decoder_table(recs):
    rows = []
    for (Np, lp), rs in agg(recs, None).items():
        f = lambda blk, k: float(np.mean([r[blk][k] for r in rs]))  # noqa: E731
        T2e, H1, T4 = f("exact", "T2"), f("first", "H"), f("exact", "T4")
        rows.append(dict(Np=Np, lp=lp, d=rs[0]["d_h"] + rs[0]["d_g"], d_h=rs[0]["d_h"], d_g=rs[0]["d_g"], n=len(rs),
                         T2_first=f("first", "T2"), T2_self=f("self", "T2"), T2_exact=T2e, T2_mf=f("mean_field", "T2"),
                         H_first=H1, H_exact=f("exact", "H"), T4=T4, T4_err=T4 - (T2e - H1),
                         walls3_first=float(np.mean([r["first"]["walls"]["0.001"] for r in rs])),
                         walls3_exact=float(np.mean([r["exact"]["walls"]["0.001"] for r in rs])),
                         atoms=float(np.mean([r["exact"]["atoms"] for r in rs])),
                         mean_lam=f("exact", "mean_lam"), mixing=float(np.mean([r["mixing_T2"] for r in rs]))))
    return rows


def resonance_table(recs):
    rows = []
    for (Np, lp), rs in agg(recs, None).items():
        for fam in ("h", "g"):
            if rs[0]["d_" + fam] == 0:
                continue
            reg = rs[0]["bookkeeping"][fam]["regime"]
            inn = [r["inner"][fam] for r in rs if r["inner"][fam].get("s2_min")]
            df = [r["dress_first"][fam] for r in rs]
            de = [r["dress_exact"][fam] for r in rs]
            rows.append(dict(Np=Np, lp=lp, fam=fam, regime=reg, d=rs[0]["d_" + fam],
                             gap=float(np.mean([x["s2_min"] / x["s2_mean"] for x in inn])) if inn else None,
                             r_low=float(np.mean([np.mean(x["r_low"][:4]) for x in inn])) if inn else None,
                             share_low1=float(np.mean([x["share_low1"] for x in inn])) if inn else None,
                             first_mean=float(np.mean([x["mean"] for x in df])), first_q50=float(np.mean([x["q50"] for x in df])),
                             first_max=float(np.max([x["max"] for x in df])),
                             exact_mean=float(np.mean([x["mean"] for x in de])), exact_max=float(np.max([x["max"] for x in de]))))
    return rows


def fig_resonance(recs3):
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6))
    ax = axs[0]
    pts = defaultdict(list)
    for r in recs3:
        for fam in ("h", "g"):
            inn = r.get("inner", {}).get(fam, {})
            if r["d_" + fam] and inn.get("s2_min"):
                pts[r["bookkeeping"][fam]["regime"]].append((r["N"], inn["s2_min"] / inn["s2_mean"]))
    for reg, mk in (("gapped", "o"), ("resonant", "s"), ("obstructed", "^")):
        if pts[reg]:
            x, y = np.array(pts[reg]).T
            ax.semilogy(x + {"gapped": -0.12, "resonant": 0.0, "obstructed": 0.12}[reg], y, mk, ms=3.5,
                        color=REG_COL[reg], label={"gapped": "gapped target", "resonant": "window-edge target",
                                                   "obstructed": "target in window"}[reg], alpha=0.85)
    Ns = np.arange(7, 16)
    ax.set_xlabel(r"old size $N$")
    ax.set_ylabel(r"$E_{\min}/\bar E$ of the inner family")
    ax.legend(fontsize=6.5, frameon=False, loc="lower left")
    ax.set_title("(a) inner gap", fontsize=8.5, loc="left")
    ax = axs[1]
    for Np, ls in ((12, ":"), (14, "--"), (16, "-")):
        for reg in ("resonant", "gapped"):
            s = []
            for r in recs3:
                if r["Np"] != Np or not r["d_h"] + r["d_g"]:
                    continue
                for fam in ("h", "g"):
                    if r["d_" + fam] and r["bookkeeping"][fam]["regime"] == reg and abs(r["lp"]) <= 2:
                        s += r["sigma2_" + fam]
            if len(s) < 20:
                continue
            s = np.sort(np.array(s))
            cc = 1 - np.arange(len(s)) / len(s)          # fraction of all states (atoms included) above K
            keep = s > 1e-8
            ax.loglog(s[keep], cc[keep], ls=ls, color=REG_COL[reg], lw=1.1, label=f"N'={Np} {reg}")
    ct = chord_tail()
    for Np, ls in ((12, ":"), (14, "--"), (16, "-")):
        c = ct.get(str(Np))
        if c:
            K, P = np.array(c["K"]), np.array(c["P"])
            m = (K >= 10) & (P > 0)
            ax.loglog(K[m], P[m], color="k", lw=0.8, ls=ls, alpha=0.7)
    ax.text(1.2e3, 0.06, r"chord: $\propto K^{-1/2}$", fontsize=7.5)
    ax.set_xlabel(r"first-order dressing cost $K=\sigma^2$")
    ax.set_ylabel(r"$P(\sigma^2 > K)$")
    ax.set_xlim(1e-2, 1e6)
    ax.set_ylim(1e-3, 1.2)
    ax.legend(fontsize=5.6, frameon=False, loc="lower left")
    ax.set_title("(b) dressing-cost tails", fontsize=8.5, loc="left")
    ax = axs[2]
    for Np, col in ((12, "#1b9e77"), (14, "#7570b3"), (16, "#d95f02")):
        rs = [r for r in recs3 if r["Np"] == Np and r["lp"] == -2]
        if not rs:
            continue
        le = np.sort(np.concatenate([r["lam"] for r in rs]))
        l1 = np.sort(np.concatenate([r["lam_first"] for r in rs]))
        cc = np.arange(1, len(le) + 1) / len(le)
        ax.semilogx(np.maximum(l1, 1e-7), cc, ls="--", color=col, lw=1.0)
        ax.semilogx(np.maximum(le, 1e-7), cc, ls="-", color=col, lw=1.2, label=f"N'={Np}")
    ax.set_xlabel(r"decoder eigenvalue $\lambda$ ($\ell'=-2$)")
    ax.set_ylabel("cumulative fraction")
    ax.set_xlim(1e-6, 1)
    ax.set_ylim(0, 0.25)
    ax.legend(fontsize=6.5, frameon=False, loc="upper left", title="exact (—), first order (- -)", title_fontsize=6)
    ax.set_title("(c) resonant sector: low tail", fontsize=8.5, loc="left")
    fig.tight_layout(pad=0.4, w_pad=0.6)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"fig_resonance.{ext}"), dpi=200)
    plt.close(fig)


def fig_decoder_sectors(recs3):
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.6), sharey=True)
    bins = np.linspace(0, 1, 41)
    for ax, Np, lps in ((axs[0], 13, (-3, -1, 1)), (axs[1], 14, (-2, 0, 2))):
        for lp, col in zip(lps, ("#c51b1b", "#6b6b6b", "#e08214") if Np == 13 else ("#e08214", "#6b6b6b", "#3b6fb6")):
            rs = [r for r in recs3 if r["Np"] == Np and r["lp"] == lp and r["d_h"] + r["d_g"]]
            if not rs:
                continue
            le = np.concatenate([r["lam"] for r in rs])
            l1 = np.concatenate([r["lam_first"] for r in rs])
            if len(le) < 10:
                for x in le:
                    ax.axvline(x, color=col, lw=1.2, ls="-", alpha=0.8)
                ax.plot([], [], color=col, label=rf"$\ell'={lp}$: exceptional state (one per seed)")
                continue
            ax.hist(le, bins=bins, density=True, histtype="step", color=col, lw=1.3, label=rf"$\ell'={lp}$ exact")
            ax.hist(l1, bins=bins, density=True, histtype="step", color=col, lw=0.8, ls="--")
        ax.set_xlabel(r"decoder eigenvalue $\lambda$")
        ax.set_title(f"N' = {Np} (dashed: first order)", fontsize=8.5, loc="left")
        ax.legend(fontsize=6.5, frameon=False, loc="upper center")
    axs[0].set_ylabel("density")
    fig.tight_layout(pad=0.4)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"fig_decoder_sectors.{ext}"), dpi=200)
    plt.close(fig)


def fig_T4(rows):
    fig, ax = plt.subplots(figsize=(3.5, 2.8))
    for r in rows:
        if r["d"] < 20:
            continue
        reg = "resonant" if abs(r["lp"]) == 2 else "gapped"
        ax.plot(r["T2_exact"] - r["H_first"], r["T4"], "o" if r["Np"] % 2 == 0 else "s", ms=3.5, color=REG_COL[reg])
    x = np.linspace(0.25, 0.85, 10)
    ax.plot(x, x, "k-", lw=0.7)
    ax.set_xlabel(r"$T_2 - H_{\rm first}$")
    ax.set_ylabel(r"$T_4$ (exact)")
    ax.plot([], [], "o", color=REG_COL["gapped"], ms=3.5, label=r"central $\ell'$")
    ax.plot([], [], "o", color=REG_COL["resonant"], ms=3.5, label=r"resonant $|\ell'|=2$")
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout(pad=0.4)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"fig_T4_sectors.{ext}"), dpi=200)
    plt.close(fig)


def main():
    os.makedirs(FIGS, exist_ok=True)
    recs3, recs5 = load(3), load(5)
    summary = {"checks_p3": checks(recs3), "checks_p5": checks(recs5)}
    print("checks q=3:", {k: v for k, v in summary["checks_p3"].items()})
    print("checks q=5:", {k: v for k, v in summary["checks_p5"].items()})
    rows = decoder_table(recs3)
    summary["decoder_p3"] = rows
    print("\n N' l'   d (h+g)   n | T2 first  self   exact  mf    | H first exact | T4     T4-(T2-H1) | walls1e-3 first exact | <lam>  mixing")
    for r in rows:
        print(f"{r['Np']:3d} {r['lp']:+d} {r['d']:5d} ({r['d_h']}+{r['d_g']}) {r['n']} | {r['T2_first']:.4f} {r['T2_self']:.4f} "
              f"{r['T2_exact']:.4f} {r['T2_mf']:.4f} | {r['H_first']:.4f} {r['H_exact']:.4f} | {r['T4']:.4f} {r['T4_err']:+.4f} | "
              f"{r['walls3_first']:.3f} {r['walls3_exact']:.3f} | {r['mean_lam']:.3f} {r['mixing']:.4f}")
    res = resonance_table(recs3)
    summary["resonance_p3"] = res
    print("\n N' l' fam regime      d    | Emin/Emean  r_low  share_low1% | first: mean  q50   max      | exact: mean  max")
    for r in res:
        g = f"{r['gap']:.1e}" if r["gap"] is not None else "   -   "
        rl = f"{r['r_low']:.3f}" if r["r_low"] is not None else "  -  "
        sh = f"{r['share_low1']:.2f}" if r["share_low1"] is not None else " - "
        print(f"{r['Np']:3d} {r['lp']:+d} {r['fam']}  {r['regime']:<10} {r['d']:5d} | {g}  {rl}  {sh} | {r['first_mean']:9.3f} "
              f"{r['first_q50']:.3f} {r['first_max']:10.3g} | {r['exact_mean']:9.3f} {r['exact_max']:10.3g}")
    ct = chord_tail()
    tails = []
    print("\n resonant tail: Hill alpha (K > 10, K > 100), P*sqrt(K) at K = 1e2: data vs chord")
    for Np in (8, 10, 12, 14, 16):
        s2 = [x for r in recs3 if r["Np"] == Np and r["lp"] == -2 and r["d_h"] for x in r["sigma2_h"]]
        s2 += [x for r in recs3 if r["Np"] == Np and r["lp"] == 2 and r["d_g"] for x in r["sigma2_g"]]
        if not s2:
            continue
        s2 = np.array(s2)
        a10, n10 = hill(s2, 10.0)
        a100, n100 = hill(s2, 100.0)
        pd = float(np.mean(s2 > 100) * 10)
        c = ct.get(str(Np))
        pc = float(np.interp(2.0, np.log10(c["K"]), c["P"]) * 10) if c else None
        tails.append(dict(Np=Np, n=len(s2), alpha_K10=a10, n_K10=n10, alpha_K100=a100, n_K100=n100,
                          Psqrt_data_K100=pd, Psqrt_chord_K100=pc))
        print(f"  N'={Np}: states {len(s2)}  alpha(K>10) = {a10:.3f} (n={n10})  alpha(K>100) = {a100:.3f} (n={n100})  "
              f"P*sqrt(K)|K=100: data {pd:.3f}  chord {pc if pc is None else round(pc, 3)}")
    summary["resonant_tail"] = tails
    summary["decoder_p5"] = decoder_table(recs5)
    with open(os.path.join(DATA, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1, default=str)
    fig_resonance(recs3)
    fig_decoder_sectors(recs3)
    fig_T4(rows)
    print("\nwrote", os.path.join(DATA, "summary.json"), "and figures in", FIGS)


if __name__ == "__main__":
    main()
