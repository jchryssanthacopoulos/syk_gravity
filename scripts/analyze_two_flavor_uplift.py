#!/usr/bin/env python3
r"""
Analysis and figures for the two-flavor pair uplift (docs/derivations.md D8; data from scripts/run_two_flavor_uplift.py).

Checks (printed and written to results/data/two_flavor_2026-10-23/summary.json):
  dims   two-step dimension formula and four-slot injectivity; E2 page vs dim B_new (d2 sectors)
  cv     the two room laws  rk D_T = min(d_new, R_T)  and  dim(B cap im D_T) = max(0, b + rk D_T - R_T);
         CV-Fig.-6-style table of U_T (fortuitous; min over seeds) with the predicted category (1 / (0,1) / 0)
  bare   bare states vs the generic count (+ tower), arithmetic walls; tower slot weights
  first  exact vs first-order bare fraction; interference statistics
  tower  kappa_00, kappa_11 of the tower vs CV's bounds, and the criterion  kappa = bound  <=>  ker D_T = 0

Figures (results/figures/two_flavor_2026-10-23/):
  fig_cv_fidelity   U_T of the fortuitous subspace by sector, channel and N, with the room-count prediction
  fig_tower_costs   tower dressing costs vs N against CV's bounds (filled: ker D_T = 0)
  fig_bare_counts   bare states vs generic count; the tower excess
  fig_first_order   first-order vs exact bare fraction (central sectors) and the interference distribution

Usage: .venv/bin/python scripts/analyze_two_flavor_uplift.py [--no-figs]
"""
import argparse
import glob
import json
import os
from collections import defaultdict

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DATA = os.path.join(ROOT, "results/data/two_flavor_2026-10-23")
FIGS = os.path.join(ROOT, "results/figures/two_flavor_2026-10-23")
CH = ("00", "10", "01", "11")
CH_LABEL = {"00": r"00 (charge-preserving)", "10": r"10 ($+\psi$)", "01": r"01 ($+\chi$)", "11": r"11 ($+\psi+\chi$)"}
# CV Fig. 6 rows: psi-heavy, chi-heavy, diagonal
CV_ROWS = [(2, 1), (3, 1), (3, 2), (1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (1, 1), (2, 2), (3, 3)]


def load(task, model):
    out = []
    for f in sorted(glob.glob(os.path.join(DATA, f"{task}_N*_{model}*.jsonl"))):
        out += [json.loads(l) for l in open(f)]
    return out


# ------------------------------------------------------------------------------------------------ dims / e2
def check_dims(model):
    recs = load("dims", model)
    e2 = {(r["N"], r["seed"], r["A"], r["B"]): r for r in load("e2", model)}
    res = dict(n=len(recs), M_mismatch=0, new_mismatch=0, inj_fail=0, Ns=sorted({r["N"] for r in recs}), d2=defaultdict(list))
    for r in recs:
        res["M_mismatch"] += r["M_pred"] != r["M"]
        res["new_mismatch"] += r["new_pred"] != r["new"]
        res["inj_fail"] += r["injectivity_deficit"] != 0
        k = (r["N"], r["seed"], r["A"], r["B"])
        if k in e2 and e2[k]["E2_total"] != r["new"]:
            res["d2"][r["N"]].append(dict(seed=r["seed"], A=r["A"], B=r["B"], E2=e2[k]["E2"], d_new=r["new"]))
    res["d2"] = {N: v for N, v in res["d2"].items()}
    return res


# ------------------------------------------------------------------------------------------------ cv
def cv_category_pred(rec, ch, symmetrized):
    """Predicted U_fort category from the counts: '1', '0' or 'mid' (generic position; D8 sec. 6)."""
    c = rec[ch]
    b, f, rk, R = rec["b_old"], rec["fort"], c["rkD"], c["R"]
    if f == 0:
        return None
    if rk >= R:                                   # D_T fills its room: B inside im D_T
        return "1"
    protected = symmetrized and rec["a"] == rec["b"] and ch in ("00", "11")
    lost = (b - rk) if protected else (f - rk)    # generic dim of F cap (im D_T)^perp
    return "0" if lost > 0 else "mid"


def cv_category_obs(U):
    if U is None:
        return None
    return "1" if U > 1 - 1e-7 else ("0" if U < 1e-7 else "mid")


def check_cv(model):
    recs = load("cv", model)
    sym = model == "sym"
    res = dict(n=0, rk_law_fail=0, dimA_law_fail=0, tower_protected_fail=0, cat_fail=0, Ns=sorted({r["N"] for r in recs}))
    table = defaultdict(lambda: defaultdict(list))
    for r in recs:
        for ch in CH:
            c = r[ch]
            res["n"] += 1
            res["rk_law_fail"] += c["rkD"] != min(c["d_new"], c["R"])
            res["dimA_law_fail"] += c["dimA"] != max(0, r["b_old"] + c["rkD"] - c["R"])
            if sym and r["a"] == r["b"] and ch in ("00", "11") and c["U_tow"] is not None:
                res["tower_protected_fail"] += c["U_tow"] < 1 - 1e-8
            table[(r["N"], r["a"], r["b"])][ch].append(r)
    # per-cell min over seeds, categories
    cells = {}
    for key, d in table.items():
        cells[key] = {}
        for ch, rs in d.items():
            Uf = [x[ch]["U_fort"] for x in rs if x[ch]["U_fort"] is not None]
            Ut = [x[ch]["U_tow"] for x in rs if x[ch]["U_tow"] is not None]
            preds = {cv_category_pred(x, ch, sym) for x in rs}
            obs = cv_category_obs(min(Uf)) if Uf else None
            # observed per-seed categories must match per-seed predictions
            for x in rs:
                po, pp = cv_category_obs(x[ch]["U_fort"]), cv_category_pred(x, ch, sym)
                if po is not None and pp is not None and po != pp:
                    res["cat_fail"] += 1
                    res.setdefault("cat_fail_cases", []).append(dict(N=x["N"], seed=x["seed"], a=x["a"], b=x["b"], ch=ch,
                                                                     obs=po, pred=pp, U=x[ch]["U_fort"]))
            cells[key][ch] = dict(U_fort_min=min(Uf) if Uf else None, U_tow_min=min(Ut) if Ut else None,
                                  pred=sorted(p for p in preds if p is not None), obs=obs, nseeds=len(rs))
    res["cells"] = {f"{k[0]}:{k[1]},{k[2]}": v for k, v in cells.items()}
    return res, cells


def print_cv_table(cells, Ns):
    for ch in CH:
        print(f"\n  channel {ch}: U_fort (min over seeds) [pred]  | columns N->N+1 = {', '.join(f'{N}->{N+1}' for N in Ns)}")
        for (a, b) in CV_ROWS:
            row = []
            for N in Ns:
                c = cells.get((N, a, b), {}).get(ch)
                if c is None or c["U_fort_min"] is None:
                    row.append("    -    ")
                else:
                    tag = {"1": "1", "0": "0", "mid": "~"}["/".join(c["pred"])] if len(c["pred"]) == 1 else "?"
                    row.append(f"{c['U_fort_min']:5.2f}[{tag}]")
            print(f"    ({a},{b}) " + " ".join(row))


# ------------------------------------------------------------------------------------------------ bare
def check_bare(model):
    recs = load("bare", model)
    sym = model == "sym"
    res = dict(n=len(recs), Ns=sorted({r["N"] for r in recs}), deviations=[], subkernel_coincidences=[],
               tower_excess=0, tower_weight_err=0.0, tower_bare_err=0.0)
    for r in recs:
        devx = r["bare"] - r["gen_exact_sub"]
        tower_excess = sym and r["A"] == r["B"] and 1 <= r["A"] <= r["N"] and r["gen_exact_sub"] == 0
        res["tower_excess"] += tower_excess and devx == 1
        if r["n01_exact"] != r["n01"] or r["n10_exact"] != r["n10"]:
            res["subkernel_coincidences"].append(dict(N=r["N"], seed=r["seed"], A=r["A"], B=r["B"], n01=r["n01"],
                                                      n01_exact=r["n01_exact"], n10=r["n10"], n10_exact=r["n10_exact"]))
        if devx != (1 if tower_excess else 0):
            res["deviations"].append(dict(N=r["N"], seed=r["seed"], A=r["A"], B=r["B"], bare=r["bare"], gen=r["gen"],
                                          gen_exact_sub=r["gen_exact_sub"], n01=r["n01"], n01_exact=r["n01_exact"],
                                          n10=r["n10"], n10_exact=r["n10_exact"]))
        if sym and "tower_slot_weights" in r and 0 <= r["A"] <= r["N"] + 1:
            N, n = r["N"], r["A"]
            w = r["tower_slot_weights"]
            res["tower_weight_err"] = max(res["tower_weight_err"], abs(w["00"] - (N + 1 - n) / (N + 1)), abs(w["11"] - n / (N + 1)))
            res["tower_bare_err"] = max(res["tower_bare_err"], abs(r["tower_bare_fraction"] - 1))
    return res, recs


# ------------------------------------------------------------------------------------------------ first
def check_first(model):
    recs = load("first", model)
    by = defaultdict(list)
    for r in recs:
        by[(r["N"], r["A"], r["B"])].append(r)
    rows = []
    for (N, A, B), rs in sorted(by.items()):
        if rs[0]["d_new"] < 8:
            continue
        rows.append(dict(N=N, A=A, B=B, d_new=rs[0]["d_new"], nseeds=len(rs),
                         mu_exact=float(np.mean([x["mu_exact_mean"] for x in rs])),
                         mu_first=float(np.mean([x["mu_first_mean"] for x in rs])),
                         H_exact=float(np.mean([x["H_exact"] for x in rs])), H_first=float(np.mean([x["H_first"] for x in rs])),
                         walls_exact=float(np.mean([x["n_one_exact"] for x in rs])), walls_first=float(np.mean([x["n_one_first"] for x in rs])),
                         interf_median=float(np.median(np.concatenate([x["interference_over_incoherent"] for x in rs])))
                         if any(len(x["interference_over_incoherent"]) for x in rs) else None))
    return rows, by


# ------------------------------------------------------------------------------------------------ tower
def check_tower():
    recs = load("tower", "sym")
    res = dict(n=len(recs), criterion_fail=0, above_bound=0)
    for r in recs:
        at = abs(r["kappa"] - r["bound"]) < 1e-9 * max(1, r["bound"])
        res["criterion_fail"] += at != (r["kerD"] == 0)
        res["above_bound"] += r["kappa"] > r["bound"] + 1e-9
    return res, recs


# ------------------------------------------------------------------------------------------------ closed forms
def check_closed(model):
    recs = load("closed", model)
    res = dict(n=len(recs), Ns=sorted({r["N"] for r in recs}), hpl_fail=0, adm_fail=0, recon=0.0, bps=0.0, dec_bare=0.0,
               dec_cv=0.0, reduced=0, leak_free=0, accretive=0, herm_bound_viol=0, sigma_min_M=np.inf, worst=None,
               refined_lt2=0)
    for r in recs:
        res["hpl_fail"] += r["hpl_dim"] != r["d_new"]
        res["adm_fail"] += r["adm_dim"] != r["d_new"]
        res["recon"] = max(res["recon"], r.get("reconstruction_err", 0.0))
        res["bps"] = max(res["bps"], r.get("bps_residual", 0.0))
        res["dec_bare"] = max(res["dec_bare"], r.get("decoder_bare_err", 0.0))
        res["dec_cv"] = max(res["dec_cv"], r.get("decoder_cv_err", 0.0))
        if "sigma_min_M" in r:
            res["reduced"] += 1
            res["leak_free"] += r["leak_norm"] < 1e-10
            res["accretive"] += r["herm_min"] > 0
            res["herm_bound_viol"] += r["herm_min"] < r["accretive_bound"] - 1e-9
            res["refined_lt2"] += r.get("refined_leak", 0.0) < 2
            if r["sigma_min_M"] < res["sigma_min_M"]:
                res["sigma_min_M"] = r["sigma_min_M"]
                res["worst"] = dict(N=r["N"], seed=r["seed"], A=r["A"], B=r["B"], mu_min=r.get("mu_min"),
                                    leak=r["leak_norm"], refined=r.get("refined_leak"))
    return res, recs


# ------------------------------------------------------------------------------------------------ figures
def fig_cv_fidelity(plt, cells, Ns):
    from matplotlib.patches import Rectangle
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 10.2))
    axes = axes.ravel()
    cmap = plt.get_cmap("RdYlGn")
    for ax, ch in zip(axes, CH):
        for i, (a, b) in enumerate(CV_ROWS):
            for j, N in enumerate(Ns):
                c = cells.get((N, a, b), {}).get(ch)
                y = len(CV_ROWS) - 1 - i
                if c is None or c["U_fort_min"] is None:
                    ax.text(j, y, "–", ha="center", va="center", fontsize=8, color="0.6")
                    continue
                U = max(0.0, c["U_fort_min"])
                ax.add_patch(Rectangle((j - 0.5, y - 0.5), 1, 1, color=cmap(U), alpha=0.85, lw=0))
                ok = len(c["pred"]) == 1 and c["pred"][0] == c["obs"]
                sym = {"1": "■", "mid": "◩", "0": "□"}.get(c["pred"][0] if len(c["pred"]) == 1 else "", "?")
                lab = f"{U:.2f}" if (U >= 0.005 or U < 1e-7) else "<.01"
                ax.text(j, y + 0.17, lab, ha="center", va="center", fontsize=8.5)
                ax.text(j, y - 0.24, sym if ok else "✗", ha="center", va="center", fontsize=8,
                        color="k" if ok else "r")
        ax.set_xlim(-0.5, len(Ns) - 0.5)
        ax.set_ylim(-0.5, len(CV_ROWS) - 0.5)
        ax.set_xticks(range(len(Ns)))
        ax.set_xticklabels([f"{N}→{N+1}" for N in Ns], fontsize=9)
        ax.set_title(CH_LABEL[ch], fontsize=11)
        ax.set_xlabel(r"$N\to N+1$", fontsize=10)
        for yy in (2.5, 7.5):
            ax.axhline(yy, color="k", lw=0.8)
        ax.set_yticks(range(len(CV_ROWS)))
        ax.set_yticklabels([f"({a},{b})" for a, b in CV_ROWS[::-1]], fontsize=9)
    for ax in axes[::2]:
        ax.set_ylabel(r"old sector $(N_\psi,N_\chi)$:  diagonal | $\chi$-heavy | $\psi$-heavy", fontsize=9.5)
    fig.suptitle(r"$U_T$ of the fortuitous subspace (CV eq. 89, min over seeds)" "\n"
                 r"prediction from the room count: ■ $U=1$ ($\mathrm{rk}\,D_T=R_T$),  ◩ $0<U<1$,  □ $U=0$;  ✗ disagreement",
                 fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig


def fig_tower_costs(plt, recs):
    seed0 = min(x["seed"] for x in recs)
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.9))
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for ax, ch in zip(axes, ("00", "11")):
        for k, n in enumerate((1, 2, 3, 4)):
            rs = sorted([r for r in recs if r["channel"] == ch and r["n"] == n and r["seed"] == seed0 and r["N"] >= n],
                        key=lambda r: r["N"])
            if not rs:
                continue
            Nx = np.array([r["N"] for r in rs])
            bound = np.array([r["bound"] for r in rs])
            ax.plot(Nx, bound, "--", color=colors[k], lw=1, alpha=0.8)
            at = np.array([r["kerD"] == 0 for r in rs])
            kap = np.array([r["kappa"] for r in rs])
            ax.plot(Nx[at], kap[at], "o", color=colors[k], ms=6.5, label=rf"$n={n}$")
            ax.plot(Nx[~at], kap[~at], "o", mfc="white", color=colors[k], ms=6.5)
        ax.set_xlabel(r"$N$  (uplift $N\to N+1$)")
        ax.set_xticks(sorted({r["N"] for r in recs}))
        ax.set_ylabel(rf"$\kappa_{{{ch}}}$ of the tower $V^n|\Omega\rangle$")
        if ch == "00":
            ax.set_ylim(-0.03, 1.45)
            ax.set_title(r"channel 00; dashed: CV bound $n/(N+1-n)$", fontsize=10)
        else:
            ax.set_yscale("symlog", linthresh=0.02, linscale=0.4)
            ax.set_ylim(-0.003, 4.0)
            ax.set_yticks([0, 0.02, 0.1, 1])
            ax.set_yticklabels(["0", "0.02", "0.1", "1"])
            ax.set_title(r"channel 11; dashed: CV bound $(N-n)/(n+1)$", fontsize=10)
        ax.legend(frameon=False, fontsize=8, loc="upper left" if ch == "11" else "upper right",
                  title=r"filled: $\ker D_T=0$, open: $\ker D_T\neq0$", title_fontsize=7.5)
    fig.tight_layout()
    return fig


def fig_bare_counts(plt, recs_sym, recs_gen):
    fig, ax = plt.subplots(figsize=(4.6, 3.7))
    for recs, mk, lab, col in ((recs_gen, "s", "generic (HTZ) couplings", "0.55"), (recs_sym, "o", "symmetrized (CCSY)", "C0")):
        x = np.array([r["gen_exact_sub"] for r in recs])
        y = np.array([r["bare"] for r in recs])
        ax.plot(x + 1, y + 1, mk, mfc="none", color=col, ms=5, label=lab, alpha=0.8)
    tw = [r for r in recs_sym if r["A"] == r["B"] and r["bare"] == r["gen_exact_sub"] + 1]
    ax.plot([r["gen_exact_sub"] + 1 for r in tw], [r["bare"] + 1 for r in tw], "*", color="C3", ms=11,
            label="excess = the tower (diagonal sectors)")
    lim = [0.8, 300]
    ax.plot(lim, lim, "k-", lw=0.6)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("generic count of bare states + 1")
    ax.set_ylabel("bare states (exact) + 1")
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    ax.set_title(r"bare states of the pair uplift, all sectors, $N\to N+1$ for $N=2,\dots,7$", fontsize=9)
    fig.tight_layout()
    return fig


def fig_first_order(plt, by_sym, Ns=(6, 7)):
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.0))
    ax = axes[0]
    styles = {6: "-", 7: ":"}
    for N in Ns:
        for k, (A, B) in enumerate(((N // 2, N // 2), (N // 2 + 1, N // 2 + 1)) if N % 2 == 0 else (((N + 1) // 2, (N + 1) // 2),)):
            rs = by_sym.get((N, A, B), [])
            if not rs:
                continue
            me = np.sort(np.concatenate([r["mu_exact"] for r in rs]))
            m1 = np.sort(np.concatenate([r["mu_first"] for r in rs]))
            c = f"C{2 * (N - 6)}"
            ax.plot(me, np.arange(1, len(me) + 1) / len(me), "-", color=c, lw=1.6,
                    label=rf"$N={N}\to{N+1}$, $({A},{B})$: exact" if k == 0 else None)
            ax.plot(m1, np.arange(1, len(m1) + 1) / len(m1), "--", color=c, lw=1.2,
                    label=rf"$N={N}\to{N+1}$, $({A},{B})$: first order" if k == 0 else None)
    ax.set_xlabel(r"bare fraction $\mu$ (eigenvalue of $G=\sum_T S_T^\dagger P_H S_T$)")
    ax.set_ylabel("cumulative fraction of BPS states")
    ax.set_title("central sector: exact vs first order (all seeds)", fontsize=10)
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    ax = axes[1]
    from matplotlib.lines import Line2D
    handles = []
    for N, col in zip(Ns, ("C0", "C2")):
        mus, vals = [], []
        for (NN, A, B), rs in by_sym.items():
            if NN == N and A == B and 1 <= A <= N:
                for r in rs:
                    me = np.array(r["mu_exact"])
                    me = me[me < 1 - 1e-8]
                    if len(me) == len(r["interference_over_incoherent"]):
                        mus += list(me)
                        vals += r["interference_over_incoherent"]
        if mus:
            ax.plot(mus, vals, ".", ms=2.2, alpha=0.35, color=col)
            handles.append(Line2D([], [], ls="", marker="o", ms=4, color=col, label=rf"$N={N}\to{N+1}$, diagonal sectors"))
    ax.plot([1.0], [-1.0], "*", color="C3", ms=13, zorder=5)
    handles.append(Line2D([], [], ls="", marker="*", ms=10, color="C3", label="the tower (exactly bare, exactly dark)"))
    ax.axhline(0, color="0.5", lw=0.6)
    ax.set_xlabel(r"exact bare fraction $\mu$ of the state")
    ax.set_ylabel(r"$(\kappa^{(1)}-\kappa_{\rm inc})/\kappa_{\rm inc}$")
    ax.set_title("interference between slots at first order", fontsize=10)
    ax.set_ylim(-1.08, 1.6)
    ax.legend(handles=handles, frameon=False, fontsize=7.5, loc="upper right")
    fig.tight_layout()
    return fig


def chi_heavy_rows(model="sym", seed=1, rows=((1, 2), (2, 3), (2, 4), (3, 4), (3, 5))):
    """LaTeX rows for the chi-heavy mechanism table: d_new, room, codim = b - dim A, U_fort for channels 00 and 10."""
    recs = {(r["N"], r["a"], r["b"]): r for r in load("cv", model) if r["seed"] == seed}
    out = []
    for N in sorted({k[0] for k in recs}):
        for (a, b) in rows:
            r = recs.get((N, a, b))
            if r is None or r["fort"] == 0:
                continue
            cells = []
            for ch in ("00", "10"):
                c = r[ch]
                U = c["U_fort"]
                Us = "1" if U > 1 - 1e-7 else ("0" if U < 1e-7 else (f"{U:.2f}" if U >= 0.005 else "$<$.01"))
                cells += [str(c["d_new"]), str(c["R"]), str(r["b_old"] - c["dimA"]), Us]
            out.append(f"{N} & $({a},{b})$ & " + " & ".join(cells[:4]) + " & " + " & ".join(cells[4:]) + r"\\")
    return out


def dims_coupling_independence():
    vals = defaultdict(set)
    for task in ("dims", "bare"):
        for model in ("sym", "generic"):
            for r in load(task, model):
                vals[(r["N"], r["A"], r["B"])].add(r.get("new", r.get("d_new")))
    return dict(sectors=len(vals), Ns=sorted({k[0] for k in vals}),
                differing=[k for k, v in vals.items() if len(v) > 1])


def hard_edge_stats():
    """Fraction of near-bare states (mu > 0.99, excluding exact walls) in diagonal sectors vs the bare-map balance."""
    first = load("first", "sym")
    bare = {(r["N"], r["seed"], r["A"], r["B"]): r for r in load("bare", "sym")}
    out = []
    by = defaultdict(list)
    for r in first:
        by[(r["N"], r["A"], r["B"])].append(r)
    for (N, A, B), rs in sorted(by.items()):
        if A != B or not (1 <= A <= N) or rs[0]["d_new"] < 20:
            continue
        mu = np.concatenate([np.array(x["mu_exact"]) for x in rs])
        near = np.mean((mu > 0.99) & (mu < 1 - 1e-9))
        b = bare.get((N, rs[0]["seed"], A, B))
        bal = None
        if b is not None:
            sl, rm = b["slots"], b["rooms"]
            unknowns = sl["00"] + sl["11"] + b["n01_exact"] + b["n10_exact"]
            bal = (unknowns, rm[4] + rm[5])
        out.append(dict(N=N, A=A, B=B, d_new=rs[0]["d_new"], nseeds=len(rs), near_bare_frac=float(near), balance=bal))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-figs", action="store_true")
    a = ap.parse_args()
    summary = {}
    for model in ("sym", "generic"):
        d = check_dims(model)
        print(f"[dims {model}] N={d['Ns']}: {d['n']} sector records; two-step mismatches M {d['M_mismatch']}, new {d['new_mismatch']}; "
              f"injectivity failures {d['inj_fail']}; d2 sectors: " + "; ".join(f"N={N}: " + ", ".join(f"({x['A']},{x['B']}) E2={x['E2']} vs {x['d_new']}" for x in v if x['seed'] == min(y['seed'] for y in v)) for N, v in d["d2"].items()))
        summary[f"dims_{model}"] = d
        cv, cells = check_cv(model)
        print(f"[cv {model}] N={cv['Ns']}: {cv['n']} (sector, channel, seed) records; rank-law failures {cv['rk_law_fail']}, "
              f"dimA-law failures {cv['dimA_law_fail']}, tower 00/11 not exact {cv['tower_protected_fail']}, "
              f"category mismatches {cv['cat_fail']}" + (": " + str(sorted({(c['N'], c['a'], c['b'], c['ch'], c['obs'], c['pred']) for c in cv.get('cat_fail_cases', [])})) if cv['cat_fail'] else ""))
        if model == "sym":
            print_cv_table(cells, cv["Ns"])
            cells_sym = cells
        summary[f"cv_{model}"] = {k: v for k, v in cv.items() if k != "cells"}
        summary[f"cv_{model}"]["cells"] = cv["cells"]
        b, brecs = check_bare(model)
        sk = sorted({(x["N"], x["A"], x["B"]) for x in b["subkernel_coincidences"]})
        print(f"\n[bare {model}] N={b['Ns']}: {b['n']} records; deviations from  bare = generic count (exact sub-kernels) + tower: "
              f"{len(b['deviations'])}; tower excess cases {b['tower_excess']}; sub-kernel coincidences (N,A,B): {sk}; "
              f"tower slot-weight error {b['tower_weight_err']:.1e}, bare-fraction error {b['tower_bare_err']:.1e}")
        for x in b["deviations"][:12]:
            print("      deviation:", x)
        summary[f"bare_{model}"] = b
        if model == "sym":
            brecs_sym = brecs
        else:
            brecs_gen = brecs
        rows, by = check_first(model)
        summary[f"first_{model}"] = rows
        if model == "sym":
            by_sym = by
            print("\n[first sym] sector-averaged bare fraction (exact / first order), H = <y(1-y)>, y = (2mu-1)^2:")
            for r in rows:
                if r["N"] >= 5 and r["d_new"] >= 50:
                    print(f"     N={r['N']} ({r['A']},{r['B']}) d={r['d_new']:4d} seeds={r['nseeds']}: mu {r['mu_exact']:.4f} / {r['mu_first']:.4f}"
                          f"  H {r['H_exact']:.4f} / {r['H_first']:.4f}  walls {r['walls_exact']:.1f} / {r['walls_first']:.1f}"
                          f"  interference median {r['interf_median']:+.3f}")
    ci = dims_coupling_independence()
    print(f"\n[dims coupling independence] {ci['sectors']} sectors, N={ci['Ns']}: differing {ci['differing']}")
    summary["dims_coupling_independence"] = ci
    print("\n[chi-heavy table rows, seed 1]")
    for l in chi_heavy_rows():
        print("   ", l)
    he = hard_edge_stats()
    print("\n[hard edge] near-bare fraction (0.99 < mu < 1) in diagonal sectors; bare-map balance (unknowns, constraints):")
    for x in he:
        print(f"    N={x['N']} ({x['A']},{x['B']}) d={x['d_new']} seeds={x['nseeds']}: near-bare {x['near_bare_frac']:.3f}  balance {x['balance']}")
    summary["hard_edge"] = he
    for model in ("sym", "generic"):
        c, _ = check_closed(model)
        if c["n"]:
            print(f"\n[closed {model}] N={c['Ns']}: {c['n']} sector records; HPL dim failures {c['hpl_fail']}, admissible-dim "
                  f"failures {c['adm_fail']}; reconstruction {c['recon']:.1e}, BPS residual {c['bps']:.1e}, decoder pencils "
                  f"{c['dec_bare']:.1e} / {c['dec_cv']:.1e}; reduced systems {c['reduced']}: leak-free {c['leak_free']}, "
                  f"Herm(M) > 0 in {c['accretive']} (refined criterion < 2 in {c['refined_lt2']}), bound violations "
                  f"{c['herm_bound_viol']}; smallest sigma_min(M) {c['sigma_min_M']:.3g} at {c['worst']}")
            summary[f"closed_{model}"] = c
    t, trecs = check_tower()
    print(f"\n[tower] {t['n']} records: criterion (kappa = bound <=> ker D_T = 0) failures {t['criterion_fail']}, above bound {t['above_bound']}")
    summary["tower"] = t
    with open(os.path.join(DATA, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1, default=lambda o: o if not isinstance(o, np.generic) else o.item())
    if a.no_figs:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.labelsize": 10, "legend.fontsize": 8, "font.family": "DejaVu Sans"})
    os.makedirs(FIGS, exist_ok=True)
    for name, fig in (("fig_cv_fidelity", fig_cv_fidelity(plt, cells_sym, summary["cv_sym"]["Ns"])),
                      ("fig_tower_costs", fig_tower_costs(plt, trecs)),
                      ("fig_bare_counts", fig_bare_counts(plt, brecs_sym, brecs_gen)),
                      ("fig_first_order", fig_first_order(plt, by_sym, Ns=tuple(N for N in (6, 7) if any(k[0] == N for k in by_sym))))):
        for ext in ("pdf", "png"):
            fig.savefig(os.path.join(FIGS, f"{name}.{ext}"), dpi=200 if ext == "png" else None)
        plt.close(fig)
        print("wrote", name)


if __name__ == "__main__":
    main()
