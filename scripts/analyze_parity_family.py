#!/usr/bin/env python3
r"""
Analyze the parity-family decoder data (scripts/run_parity_family.py).  Prints markdown tables and writes
results/data/parity_family_2026-09-29/summary.csv and figures in results/figures/parity_family_2026-09-29/.

Per (N', model, s): mean +- s.e. over realizations of
  r = rho_2, effective parameters a_n = rho_n^{1/(n-1)} (free compression: a_n = a), ratios rho_n / a^{n-1},
  mu-moments T_m (m = 2, 4, 6) vs free compression Wachter(a, b_S), T_4 excess over free, KS to Wachter, atoms,
  and the double-scaled chord prediction for T_2 (sum_n P_n x^{2n}; BLY) with x exact and x double-scaled.
"""
import argparse
import csv
import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from dssyk_n2 import hh_length_distribution  # noqa: E402
from nullmodels import free_compression_mu_moments  # noqa: E402


def ms(v):
    v = np.asarray(v, float)
    return v.mean(), (v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else float("nan"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=os.path.join(ROOT, "results/data/parity_family_2026-09-29"))
    ap.add_argument("--figs", default=os.path.join(ROOT, "results/figures/parity_family_2026-09-29"))
    a = ap.parse_args()
    rows = []
    for f in sorted(glob.glob(os.path.join(a.data, "raw_q*_N*.jsonl"))):
        rows += [json.loads(l) for l in open(f) if l.strip()]
    G = defaultdict(list)
    for r in rows:
        for e in r["family"]:
            G[(r["q"], r["Np"], r["model"], e["s"])].append((r, e))

    out = []
    print("| N′ | model | s | x_exact | x_ds | r | r/a | a₄/a | a₆/a | T₄/T₄^free | T₆/T₆^free | KS_W | atoms | T₂ chord (x_ds / x_ex) |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for key in sorted(G):
        q, N, model, s = key
        L = G[key]
        aa = L[0][0]["a"]
        e0 = L[0][1]
        r_m, r_se = ms([e["r"] for _, e in L])
        a4 = ms([(e["rho"][3]) ** (1 / 3) / aa if e["rho"][3] > 0 else np.nan for _, e in L])
        a6 = ms([(e["rho"][5]) ** (1 / 5) / aa if e["rho"][5] > 0 else np.nan for _, e in L])
        T4r, T6r = [], []
        for _, e in L:
            fm = free_compression_mu_moments(aa, e["b"], 6)
            T4r.append(e["mu_moments"][3] / fm[3])
            T6r.append(e["mu_moments"][5] / fm[5])
        T4 = ms(T4r)
        T6 = ms(T6r)
        ks = ms([e["ks_wachter"] for _, e in L])
        atoms = np.mean([(e["n0"] + e["n1"]) / L[0][0]["d"] for _, e in L])
        j = (N // 2 - N / 2) / q
        P = hh_length_distribution(2 * q * q / N, j)
        n = np.arange(len(P))
        t2_ds = float((P * e0["x_ds"] ** (2 * n)).sum())
        t2_ex = float((P * e0["x_exact"] ** (2 * n)).sum())
        o = dict(q=q, Np=N, model=model, s=s, nreal=len(L), a=aa, x_exact=e0["x_exact"], x_ds=e0["x_ds"],
                 r=r_m, r_se=r_se, a4_over_a=a4[0], a4_se=a4[1], a6_over_a=a6[0], a6_se=a6[1],
                 T4_over_free=T4[0], T4_se=T4[1], T6_over_free=T6[0], T6_se=T6[1], ks=ks[0], atoms=atoms,
                 T2_chord_xds=t2_ds, T2_chord_xex=t2_ex)
        out.append(o)
        print(f"| {N} | {model} | {s} | {e0['x_exact']:+.3f} | {e0['x_ds']:.3f} | {r_m:.4f}±{r_se:.4f} | {r_m / aa:.3f} | "
              f"{a4[0]:.3f}±{a4[1]:.3f} | {a6[0]:.3f}±{a6[1]:.3f} | {T4[0]:.4f}±{T4[1]:.4f} | {T6[0]:.4f}±{T6[1]:.4f} | "
              f"{ks[0]:.3f} | {atoms:.3f} | {t2_ds:.3f} / {t2_ex:.3f} |")
    with open(os.path.join(a.data, "summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)


if __name__ == "__main__":
    main()
