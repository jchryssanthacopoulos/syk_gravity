#!/usr/bin/env python3
r"""
Heterogeneity prediction of the decoder fourth moment from the first-order uplift (docs/derivations.md D6 §5).

Exact identity (y = mu^2 = cos^2 of twice the principal angle, T_2m = <y^m>):
    T_4 = T_2 - H,      H = <y (1 - y)>   (the "heterogeneity": weight at intermediate angles).
At equal T_2 the localization (decoder - OU, D5) is therefore  L = H_OU(x*) - H_dec.
Prediction tested here: H_dec ~= H_first, computed from the OLD model (N'-1 modes) alone through the first-order
uplift spectrum  y = ((1 - sigma^2)/(1 + sigma^2))^2  (sigma: singular values of K = Q0^+ eta^ on B_o^F, and of K' on
B_o^{F-1}).  Then
    T_4^pred = T_2 - H_first,       L^pred = T_4^pred - T_4^OU(x*).
Inputs:
  results/data/uplift_2026-10-12/walls_N*.jsonl      (scripts/run_uplift_walls.py: first order per seed, exact if run)
  results/data/ou_fourarc_2026-10-05/ou_p3_N*.jsonl  (scripts/run_ou_fourarc.py: decoder all-site T_2, T_4, OU curve)
  results/data/dssyk_2026-09-28/comparison.csv       (D2: chord prediction of r, hence T_2 with no free parameter)

Usage: .venv/bin/python scripts/analyze_uplift_heterogeneity.py [--Np 10 11 12 13 14 15 16] [--json out] [--fig out]
"""
import argparse
import csv
import json
import os
import sys

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from analyze_ou_fourarc import analyze  # noqa: E402

WALLS = os.path.join(ROOT, "results/data/uplift_2026-10-12")
CHORD = os.path.join(ROOT, "results/data/dssyk_2026-09-28/comparison.csv")


def het(lam):
    y = (2.0 * np.asarray(lam) - 1.0) ** 2
    return float(np.mean(y * (1.0 - y))), float(np.mean(y)), float(np.mean(y * y))


def load_walls(N):
    path = os.path.join(WALLS, f"walls_N{N}.jsonl")
    if not os.path.exists(path):
        return None
    rows = {}
    for line in open(path):
        r = json.loads(line)
        rows[r["seed"]] = {**rows.get(r["seed"], {}), **r}          # later records (e.g. exact) complete earlier ones
    out = dict(seeds=sorted(rows), H_first=[], H_exact=[], H_pair=[], T2_exact=[], T4_exact=[], fold_ks=[])
    for s in out["seeds"]:
        r = rows[s]
        h1 = het(r["first_lam"])[0]
        out["H_first"].append(h1)
        if "lam" in r:
            he, t2, t4 = het(r["lam"])
            out["H_exact"].append(he)
            out["H_pair"].append((h1, he))
            out["T2_exact"].append(t2)
            out["T4_exact"].append(t4)
    return out


def chord_T2():
    """T_2 = 1 - 4 b(1-b)(1 - r) with b = 1 - P/N' (D1 Step 1; exact given r), r from BLY (D2), main sector only."""
    out = {}
    for row in csv.DictReader(open(CHORD)):
        N, P, j = int(row["Np"]), int(row["P"]), float(row["j"])
        if P != N // 2:
            continue
        b = 1.0 - P / N
        out[N] = 1.0 - 4.0 * b * (1.0 - b) * (1.0 - float(row["r_chord"]))
    return out


def mean_se(v):
    v = np.asarray(v, float)
    if len(v) == 0:
        return np.nan, np.nan
    return float(v.mean()), float(v.std(ddof=1) / np.sqrt(len(v))) if len(v) > 1 else np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Np", type=int, nargs="+", default=[10, 11, 12, 13, 14, 15, 16])
    ap.add_argument("--json", default=os.path.join(WALLS, "heterogeneity_analysis.json"))
    ap.add_argument("--fig", default=os.path.join(ROOT, "results/figures/uplift_2026-10-12/fig_heterogeneity.pdf"))
    ap.add_argument("--nboot", type=int, default=300)
    ap.add_argument("--sym", default="H", help="symbol for <y(1-y)> in figure labels (the 2026-10-13 report uses W)")
    a = ap.parse_args()
    t2c = chord_T2()
    res = []
    print("H = <y(1-y)>.  'first' = first-order uplift from the old model (mean over seeds +- SE); 'exact pair' = exact "
          "decoder on the SAME realizations (last site).")
    print(f"{'N':>3} {'nseed':>5} {'H_first':>15} {'H_exact(pair)':>15} {'exact-first':>16} {'H_dec(all-site)':>16}")
    for N in a.Np:
        w = load_walls(N)
        if w is None:
            continue
        H1, H1se = mean_se(w["H_first"])
        d = [he - h1 for h1, he in w["H_pair"]]
        He, Hese = mean_se([he for _, he in w["H_pair"]])
        dm, dse = mean_se(d)
        e = dict(N=N, n_first=len(w["H_first"]), H_first=H1, H_first_se=H1se, n_exact=len(w["H_pair"]), H_exact=He,
                 H_exact_se=Hese, exact_minus_first=dm, exact_minus_first_se=dse, T2_chord=t2c.get(N, np.nan),
                 T2_exact_walls=mean_se(w["T2_exact"])[0], T4_exact_walls=mean_se(w["T4_exact"])[0])
        try:
            o = analyze(N, nboot=a.nboot)["decoder"][1]
            e.update(T2_dec=o["T2_dec"], T2_dec_se=o["T2_dec_se"], T4_dec=o["T4_dec"], T4_dec_se=o["T4_dec_se"],
                     T4_ou=o["T4_ou"], T4_ou_se=o["T4_ou_se"], T4_lib=o["T4_lib"], x_star=o["x_star"],
                     L_meas=o["split_loc"], L_meas_se=o["dT4_se"])
            e["H_dec"] = e["T2_dec"] - e["T4_dec"]
            e["T4_pred"] = e["T2_dec"] - H1
            e["T4_pred_se"] = float(np.hypot(H1se if np.isfinite(H1se) else 0.0, e["T2_dec_se"]))
            e["T4_pred_chord"] = e["T2_chord"] - H1
            e["L_pred"] = e["T4_pred"] - e["T4_ou"]
        except FileNotFoundError:
            pass
        res.append(e)
        print(f"{N:>3} {e['n_first']:>5} {H1:.4f}+-{H1se:.4f} {He:.4f}+-{Hese:.4f}  {dm:+.4f}+-{dse:.4f}  "
              f"{e.get('H_dec', np.nan):.4f}")
    print("\nT_4 of the decoder (all-site mean, OU campaign) vs predictions at the decoder's own T_2:")
    print(f"{'N':>3} {'T2_dec':>8} {'T4_dec':>16} {'T2-H_first':>16} {'err':>8} {'T2chord-H_first':>16} {'err':>8} "
          f"{'lib':>8} {'err':>8} {'OU(x*)':>8} {'err':>8} | {'L_meas':>15} {'L_pred':>8}")
    for e in res:
        if "T4_dec" not in e:
            continue
        print(f"{e['N']:>3} {e['T2_dec']:.4f} {e['T4_dec']:.4f}+-{e['T4_dec_se']:.4f} {e['T4_pred']:.4f}+-{e['T4_pred_se']:.4f} "
              f"{e['T4_pred'] - e['T4_dec']:+.4f} {e['T4_pred_chord']:16.4f} {e['T4_pred_chord'] - e['T4_dec']:+.4f} "
              f"{e['T4_lib']:.4f} {e['T4_lib'] - e['T4_dec']:+.4f} {e['T4_ou']:.4f} {e['T4_ou'] - e['T4_dec']:+.4f} | "
              f"{e['L_meas']:.4f}+-{e['L_meas_se']:.4f} {e['L_pred']:.4f}")
    os.makedirs(os.path.dirname(a.json), exist_ok=True)
    json.dump(res, open(a.json, "w"), indent=1)
    print("wrote", a.json)
    if a.fig:
        plot(res, a.fig, a.sym)


def plot(res, path, sym="H"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.labelsize": 10, "legend.fontsize": 8.2})
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.8))
    E = [e for e in res if "T4_dec" in e]
    NE = np.array([e["N"] for e in E])
    # (a) heterogeneity H = T2 - T4
    ax[0].plot(NE, [e["T2_dec"] - e["T4_lib"] for e in E], "^:", color="tab:green", mfc="none",
               label="free liberation at equal $T_2$")
    ax[0].plot(NE, [e["T2_dec"] - e["T4_ou"] for e in E], "s--", color="tab:gray", mfc="none",
               label="isotropic OU at equal $T_2$")
    ax[0].errorbar(NE, [e["H_dec"] for e in E], yerr=[e["T4_dec_se"] for e in E], fmt="ko", ms=5,
                   label="decoder (exact numerics)")
    W = [e for e in res if np.isfinite(e.get("H_exact", np.nan)) and "T4_dec" not in e]
    if W:
        ax[0].errorbar([e["N"] for e in W], [e["H_exact"] for e in W], yerr=[e["H_exact_se"] for e in W], fmt="k*",
                       ms=9, label="decoder, out of sample")
    ax[0].errorbar([e["N"] for e in res], [e["H_first"] for e in res], yerr=[e["H_first_se"] for e in res], fmt="D-",
                   color="tab:red", mfc="none", ms=5, label="first-order uplift (old model only)")
    ax[0].set_xlabel(r"$N'$")
    ax[0].set_ylabel(rf"${sym}=\langle y(1-y)\rangle = T_2-T_4$")
    ax[0].set_title("heterogeneity of decoder angles", fontsize=10)
    ax[0].legend(frameon=False, loc="lower right")
    # (b) error of T4 predictors
    def err(key):
        return [e[key] - e["T4_dec"] for e in E]
    ax[1].axhline(0, color="k", lw=0.6)
    ax[1].plot(NE, err("T4_lib"), "^:", color="tab:green", mfc="none", label="free liberation")
    ax[1].plot(NE, err("T4_ou"), "s--", color="tab:gray", mfc="none", label="isotropic OU")
    ax[1].errorbar(NE, err("T4_pred"), yerr=[e["T4_pred_se"] for e in E], fmt="D-", color="tab:red", mfc="none", ms=5,
                   label=rf"$T_2-{sym}_{{\rm first}}$")
    ax[1].plot(NE, err("T4_pred_chord"), "v-.", color="tab:blue", mfc="none",
               label=rf"$T_2^{{\rm chord}}-{sym}_{{\rm first}}$ (no $N'$-mode input)")
    if W:   # out of sample (no OU campaign): T2 - H_first vs exact T4 on the same realizations = H_exact - H_first
        ax[1].errorbar([e["N"] for e in W], [e["exact_minus_first"] for e in W],
                       yerr=[e["exact_minus_first_se"] for e in W], fmt="*", color="tab:red", ms=11,
                       label=rf"$T_2-{sym}_{{\rm first}}$, out of sample (pre-registered)")
    ax[1].set_xlabel(r"$N'$")
    ax[1].set_ylabel(r"$T_4^{\rm predicted}-T_4^{\rm decoder}$")
    ax[1].set_title("decoder fourth moment: prediction error", fontsize=10)
    ax[1].set_ylim(-0.115, 0.035)
    ax[1].set_yticks([-0.06, -0.04, -0.02, 0.0, 0.02])
    ax[1].legend(frameon=False, loc="lower center", ncol=2, fontsize=7.6, columnspacing=1.0)
    fig.tight_layout()
    fig.savefig(path)
    print("wrote", path)


if __name__ == "__main__":
    main()
