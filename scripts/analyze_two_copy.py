#!/usr/bin/env python3
r"""
Exact two-copy (fourth-moment) structure of BPS projectors (research/notes/decoder_moments_2026-09-29.md, §10).

From the exact size weights w_k of single realizations (run_size_decomposition.py) compute the U(N)-twirled
transmissions phi_k of X -> P X P via the universal realignment matrix (size_decomposition.transmission_matrix), and:
  checks:     phi_0 = a;  sum_k phi_k dim V_k = d^2;  (half filling) phi_1 / a = sum_k w_k chi_hat_k(U)  (= E T_2)
  prediction: E tau(P1 P2 P1 P2) = sum_k phi_k w_k  (two independent models)
            = (2a^2 - a^3) + sum_{k>=1} (phi_k - a^2) w_k
  measured:   two-model E[nu^2] (run_two_model_overlap.py), plain mean and control-variate estimate using the exact
              E[nu] = a (regression on nu - a).
Prints markdown; writes results/data/parity_family_2026-09-29/two_copy_tables.md.  Light (< 0.5 GB).
"""
import json
import os
import sys

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from size_decomposition import chi_hat, dim_Vk, transmission_matrix, transmissions_from_weights  # noqa: E402

DATA = os.path.join(ROOT, "results/data/parity_family_2026-09-29")


def load(name):
    f = os.path.join(DATA, name)
    return [json.loads(l) for l in open(f)] if os.path.exists(f) else []


def main():
    out = []
    p = out.append
    p("### Exact checks and transmissions (SYK, mean over realizations)\n")
    p("| N′ | n | max|φ₀−a| | max|Σφ_k dimV_k − d²|/d² | max|φ₁/a − Σw_kχ̂_k| | φ_k − a² (k = 1, 2, 3, 4, …) |")
    p("|---|---|---|---|---|---|")
    pred_rows = []
    for N in (8, 10, 12, 13, 14):
        P = N // 2
        rows = load(f"size_q3_N{N}.jsonl")
        if not rows:
            continue
        A = transmission_matrix(N, P)
        for model in ("syk", "haar"):
            R = [r for r in rows if r["model"] == model]
            if not R:
                continue
            a, d = R[0]["a"], R[0]["d"]
            phis = [transmissions_from_weights(np.array(r["w"]), d, N, P, A) for r in R]
            e0 = max(abs(ph[0] - a) for ph in phis)
            esum = max(abs(sum(ph[k] * dim_Vk(N, k) for k in range(len(ph))) - d * d) / d ** 2 for ph in phis)
            et2 = (max(abs(ph[1] / a - sum(r["w"][k] * chi_hat(N, 1, k) for k in range(len(r["w"]))))
                       for ph, r in zip(phis, R)) if N % 2 == 0 else float("nan"))
            phi = np.mean(phis, axis=0)
            w = np.mean([r["w"] for r in R], axis=0)
            if model == "syk":
                p(f"| {N} | {len(R)} | {e0:.1e} | {esum:.1e} | {et2:.1e} | " + ", ".join(f"{x:+.4f}" for x in phi[1:] - a * a) + " |")
            pred_rows.append((N, model, len(R), a, float(phi @ w), float(((phi[1:] - a * a) * w[1:]).sum())))
    p("\n### Two independent models: E τ(P₁P₂P₁P₂) predicted from single-model data vs measured\n")
    p("| N′ | model | size realizations | predicted | Σ_{k≥1}(φ_k−a²)w_k | measured: plain (pairs) | measured: control variate | free 2a²−a³ | pred/free | meas/free |")
    p("|---|---|---|---|---|---|---|---|---|---|")
    for N, model, nsz, a, pred, exc in pred_rows:
        pairs = [q for q in load(f"two_model_q3_N{N}.jsonl") if q["model"] == model]
        free = 2 * a * a - a ** 3
        if pairs:
            x = np.array([q["nu_moments"][0] for q in pairs])
            y = np.array([q["nu_moments"][1] for q in pairs])
            n = len(y)
            plain, se = y.mean(), (y.std(ddof=1) / np.sqrt(n) if n > 1 else float("nan"))
            if n > 2:
                beta = np.cov(x, y, ddof=1)[0, 1] / np.var(x, ddof=1)
                resid = y - beta * (x - a)
                cv, cse = resid.mean(), resid.std(ddof=1) / np.sqrt(n)
                cvs = f"{cv:.5f} ± {cse:.5f}"
                mf = f"{cv / free:.4f}"
            else:
                cvs, mf = "–", f"{plain / free:.4f}"
            ms = f"{plain:.5f} ± {se:.5f} ({n})"
        else:
            ms, cvs, mf = "–", "–", "–"
        p(f"| {N} | {model} | {nsz} | {pred:.5f} | {exc:+.5f} | {ms} | {cvs} | {free:.5f} | {pred / free:.4f} | {mf} |")
    txt = "\n".join(out)
    print(txt)
    with open(os.path.join(DATA, "two_copy_tables.md"), "w") as f:
        f.write(txt + "\n")


if __name__ == "__main__":
    main()
