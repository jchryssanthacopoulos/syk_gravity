#!/usr/bin/env python3
r"""
Analyze decoder T_4 data (scripts/run_decoder_T4.py; research note section 11).

Exact decomposition (per realization, any probe U in U(N)), with P' = U P U = a 1 + Y:
    T_4 = a^2 + 2a (T_2 - a) + Q_4,     Q_4 = Tr(P Y P Y)/d = S + R,   S = sum_{k>=1} phi_k w_k = c_0 - a^2
Free (Haar): T_2 = a, Q_4 = S = a^2 (1 - a)  =>  T_4 = 2a^2 - a^3.
Coherent-mixture ansatz (Y ~ eta P_{>=1} + sqrt(1-eta^2) Y_indep, eta = (T_2 - a)/(1 - a)):
    R ~ eta^2 R(1),   R(1) = (1-a)^2 - S     <=>    T_4 ~ T_2^2 + (1 - eta^2) S
Closed prediction from the size spectrum alone: T_2 -> sum_k w_k chi_hat_k(U), S from w (transmission matrix).

Prints markdown tables; writes results/data/parity_family_2026-09-29/decoderT4_tables.md.
"""
import json
import os

import sys

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from nullmodels import free_compression_mu_moments  # noqa: E402

DATA = os.path.join(ROOT, "results/data/parity_family_2026-09-29")


def main():
    out = []
    p = out.append
    p("### (A) Single-site decoder: exact decomposition of T_4 − T_4^free (mean over realizations × sites)\n\nFree baseline = Wachter(a, b) moments of mu (b = 1 − P/N′; at even N′ T₂^free = a, T₄^free = 2a² − a³, Q₄^free = a²(1−a)).\n")
    p("| N′ | model | n (real × sites) | a | T₂ | T₄ | T₄/free | 2a(T₂−T₂^free) | S − Q₄^free | R | shares of excess: linear / two-model / R |")
    p("|---|---|---|---|---|---|---|---|---|---|---|")
    rows_B = []
    for N in (8, 10, 12, 13, 14):
        f = os.path.join(DATA, f"decoderT4_q3_N{N}.jsonl")
        if not os.path.exists(f):
            continue
        rows = [json.loads(l) for l in open(f)]
        for model in ("syk", "haar"):
            R_ = [r for r in rows if r["model"] == model]
            if not R_:
                continue
            a = R_[0]["a"]
            b = 1 - (N // 2) / N
            fm = free_compression_mu_moments(a, b, 4)        # free (Wachter(a,b)) E[mu^2], E[mu^4]; b = 1/2 at even N'
            T2free, free = fm[1], fm[3]
            Q4free = free - a * a - 2 * a * (T2free - a)
            T2, T4, lin, twm, rem, mod_err, etaQ, eta2, predw = [], [], [], [], [], [], [], [], []
            for r in R_:
                S = r["c0"] - a * a
                R1 = r["R_identity"]
                for e in r["probes"]:
                    if e["kind"] != "site":
                        continue
                    t2, t4 = e["T2"], e["T4"]
                    eta = (t2 - a) / (1 - a)
                    T2.append(t2); T4.append(t4)
                    lin.append(2 * a * (t2 - T2free)); twm.append(S - Q4free); rem.append(e["R_exact"])
                    mod_err.append((t2 * t2 + (1 - eta * eta) * S) / t4 - 1)
                    etaQ.append(e["R_exact"] / R1); eta2.append(eta * eta)
                    tp = e["T2pred"]
                    ep = (tp - a) / (1 - a)
                    predw.append(tp * tp + (1 - ep * ep) * S)
            T2m, T4m = np.mean(T2), np.mean(T4)
            L, Tw, Rm = np.mean(lin), np.mean(twm), np.mean(rem)
            exc = T4m - free
            sh = f"{L / exc:.2f} / {Tw / exc:.2f} / {Rm / exc:.2f}" if abs(exc) > 1e-3 else "–"
            p(f"| {N} | {model} | {len(R_)}×{len(T2) // len(R_)} | {a:.4f} | {T2m:.4f} | {T4m:.5f} | {T4m / free:.3f} | "
              f"{L:+.5f} | {Tw:+.5f} | {Rm:+.5f} | {sh} |")
            rows_B.append((N, model, a, T4m, np.std(T4, ddof=1) / np.sqrt(len(T4)), np.mean(mod_err), np.max(np.abs(mod_err)),
                           np.mean(etaQ), np.mean(eta2), np.mean(predw), R_[0]["R_identity_check"]))
    p("\n### (B) Coherent-mixture law T₄ ≈ T₂² + (1 − η²)S, η = (T₂ − a)/(1 − a)\n")
    p("| N′ | model | R(1) = (1−a)²−S | ⟨R/R(1)⟩ | ⟨η²⟩ | model error: mean / max |·| (per realization × site) | T₄ measured (mean ± se) | T₄ from w alone |")
    p("|---|---|---|---|---|---|---|---|")
    for N, model, a, T4m, T4se, me, mx, eq, e2, pw, R1 in rows_B:
        p(f"| {N} | {model} | {R1:+.4f} | {eq:+.3f} | {e2:.3f} | {me * 100:+.2f}% / {mx * 100:.2f}% | {T4m:.5f} ± {T4se:.5f} | {pw:.5f} |")
    p("\n### (C) Parity family (s ≥ 2): model error per realization\n")
    p("| N′ | s | ⟨T₄⟩ | model error mean / max |·| |")
    p("|---|---|---|---|")
    for N in (12, 13, 14):
        f = os.path.join(DATA, f"decoderT4_q3_N{N}.jsonl")
        if not os.path.exists(f):
            continue
        rows = [json.loads(l) for l in open(f) if json.loads(l)["model"] == "syk"]
        smax = max(e["s"] for e in rows[0]["probes"])
        for s in range(2, smax + 1):
            errs, t4s = [], []
            for r in rows:
                a = r["a"]
                S = r["c0"] - a * a
                for e in r["probes"]:
                    if e["kind"] == "parity" and e["s"] == s:
                        eta = (e["T2"] - a) / (1 - a)
                        errs.append((e["T2"] ** 2 + (1 - eta * eta) * S) / e["T4"] - 1)
                        t4s.append(e["T4"])
            p(f"| {N} | {s} | {np.mean(t4s):.5f} | {np.mean(errs) * 100:+.2f}% / {np.max(np.abs(errs)) * 100:.2f}% |")
    txt = "\n".join(out)
    print(txt)
    with open(os.path.join(DATA, "decoderT4_tables.md"), "w") as fh:
        fh.write(txt + "\n")


if __name__ == "__main__":
    main()
