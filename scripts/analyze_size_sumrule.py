#!/usr/bin/env python3
r"""
Tables for the exact size sum rule (research/notes/decoder_moments_2026-09-29.md):
  (1) operator-size weights w_k of the BPS projector (SYK vs Haar), mean over realizations;
  (2) E[T_2] = sum_k w_k chi_hat_k(U) for the single-site decoder vs the campaign's mode-averaged r
      (results/data/chord_invariants_2026-09-28/summary.csv);
  (3) w_{2n} vs the double-scaled wormhole-length distribution P_n (BLY), and chi_hat_{2n} vs x^{2n};
  (4) two independent models vs free compression (two_model_q3_N*.jsonl).
Light (< 0.1 GB).
"""
import csv
import glob
import json
import os
import sys
from math import exp

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from dssyk_n2 import hh_length_distribution  # noqa: E402
from size_decomposition import chi_hat  # noqa: E402

DATA = os.path.join(ROOT, "results/data/parity_family_2026-09-29")


def main():
    camp = {int(r["Np"]): (float(r["r"]), float(r["r_se"])) for r in
            csv.DictReader(open(os.path.join(ROOT, "results/data/chord_invariants_2026-09-28/summary.csv")))
            if r["q"] == "3" and r["model"] == "syk"}
    print("### (1)-(2) size weights and the T_2 sum rule (single-site decoder)\n")
    print("| N′ | model | n | w_0 (= a) | w_1 | w_2 | w_3 | w_4 | w_5 | w_6 | Σ w_k χ̂_k(U) | campaign r |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for f in sorted(glob.glob(os.path.join(DATA, "size_q3_N*.jsonl")), key=lambda p: int(p.split("N")[-1].split(".")[0])):
        rows = [json.loads(l) for l in open(f)]
        for model in ("syk", "haar"):
            R = [r for r in rows if r["model"] == model]
            if not R:
                continue
            N = R[0]["Np"]
            w = np.mean([r["w"] for r in R], axis=0)
            pred = np.mean([r["family"][0]["T2_pred"] for r in R])
            cells = " | ".join(f"{w[k]:.4f}" if k < len(w) else "–" for k in range(7))
            c = f"{camp[N][0]:.4f}±{camp[N][1]:.4f}" if (model == "syk" and N in camp) else "–"
            print(f"| {N} | {model} | {len(R)} | {cells} | {pred:.4f} | {c} |")
    print("\n### (3) exact size weights vs BLY wormhole-length distribution (even N′, j = 0)\n")
    print("| N′ | n | w_{2n} (exact) | P_n (BLY, λ=18/N′) | χ̂_{2n}(U) (exact) | x^{2n} (x = e^{-2p/N′}) |")
    print("|---|---|---|---|---|---|")
    for N in (8, 10, 12, 14):
        f = os.path.join(DATA, f"size_q3_N{N}.jsonl")
        if not os.path.exists(f):
            continue
        w = np.mean([json.loads(l)["w"] for l in open(f) if json.loads(l)["model"] == "syk"], axis=0)
        P = hh_length_distribution(18 / N, 0.0)
        for n in range(len(w[0::2])):
            print(f"| {N} | {n} | {w[2 * n]:.4f} | {P[n]:.4f} | {chi_hat(N, 1, 2 * n):.4f} | {exp(-6 / N) ** (2 * n):.4f} |")
    print("\n### (4) two independent SYK models: nu = cos^2(principal angles) between B_1 and B_2\n")
    print("| N′ | model | pairs | E[ν]/a | E[ν²]/free | E[ν³]/free | KS to Wachter(a,a) |")
    print("|---|---|---|---|---|---|---|")
    for f in sorted(glob.glob(os.path.join(DATA, "two_model_q3_N*.jsonl")), key=lambda p: int(p.split("N")[-1].split(".")[0])):
        rows = [json.loads(l) for l in open(f)]
        for model in ("syk", "haar"):
            R = [r for r in rows if r["model"] == model]
            if not R:
                continue
            a = R[0]["a"]
            fr = [a, 2 * a * a - a ** 3, 5 * a ** 3 - 6 * a ** 4 + 2 * a ** 5]
            m = np.array([r["nu_moments"][:3] for r in R]) / np.array(fr)
            se = m.std(0, ddof=1) / np.sqrt(len(R)) if len(R) > 1 else np.full(3, np.nan)
            ks = np.mean([r["ks_wachter_aa"] for r in R])
            print(f"| {R[0]['Np']} | {model} | {len(R)} | {m[:, 0].mean():.4f}±{se[0]:.4f} | {m[:, 1].mean():.4f}±{se[1]:.4f} | "
                  f"{m[:, 2].mean():.4f}±{se[2]:.4f} | {ks:.3f} |")


if __name__ == "__main__":
    main()
