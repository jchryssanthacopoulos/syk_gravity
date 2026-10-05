#!/usr/bin/env python3
r"""
Analysis of scripts/run_ou_fourarc.py output (docs/derivations.md D5; Phase A1/A2).

A1 (OU curve).  T_2m(x) = E tau((P P_x)^m) at x = cos(L pi / K), L = 0..K/2, from sample means (SE over samples).
  - chaos coefficients: T_2(x) = sum_n pi_n y^n, T_4(x) = sum_n theta_n y^n, y = x^2, by weighted least squares
    with T(1) = 1 imposed, degree --deg (bootstrap over samples for errors);
  - remainder R^OU(x) = T_4(x) - c0 - 2a (T_2(x) - a) with c0 = T_4(0) (two-model value; S is x-independent to
    < 1e-6, checked from per-realization size data).
A2 (decoder vs OU at matched T_2).  For each probe size s = |A| (s = 1: the decoder), the OU curve is parametrized by
  T_2 (monotone in x) and T_{2m}^OU is interpolated (PCHIP in T_2) at T_2^dec.  Reports x_* (matched), x_exact(s)
  (Krawtchouk), T_4, T_6, T_8 decoder vs OU(x_*), and R_dec vs R^OU(x_*).  Errors: bootstrap over samples (the OU
  curve and the decoder resampled jointly, since they share realizations).

Usage: .venv/bin/python scripts/analyze_ou_fourarc.py --Np 8 9 10 11 12 13 14 [--json out.json]
"""
import argparse
import json
import os
import sys
from math import comb

import numpy as np
from scipy.interpolate import PchipInterpolator

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DATA = os.path.join(ROOT, "results/data/ou_fourarc_2026-10-05")


def load(N, p=3, F=None, data=DATA, site=None):
    F = N // 2 if F is None else F
    path = os.path.join(data, f"ou_p{p}_N{N}_F{F}" + ("" if site is None else f"_site{site}") + ".jsonl")
    return [json.loads(l) for l in open(path)]


def krawtchouk_x(N, p, s):
    """E_J (-1)^{|A cap J|} for a uniformly random p-subset J and |A| = s."""
    return sum((-1) ** j * comb(s, j) * comb(N - s, p - j) for j in range(0, min(s, p) + 1)) / comb(N, p)


def ou_table(rows, kmax):
    """Per-sample OU moments: returns xs (L = 0..K/2, x descending from 1), T[sample, L, k] (k = 1..kmax)."""
    K = rows[0]["K"]
    xs = np.array([1.0] + [o["x"] for o in rows[0]["ou"]])
    T = np.ones((len(rows), K // 2 + 1, kmax))
    for i, r in enumerate(rows):
        for o in r["ou"]:
            T[i, o["lag"], :] = o["moments"][:kmax]
    return xs, T


def probe_table(rows, kmax):
    """Per-sample, per-s decoder moments averaged over probes: dict s -> array[sample, k]."""
    out = {}
    for i, r in enumerate(rows):
        acc = {}
        for pr in r["probes"]:
            acc.setdefault(pr["s"], []).append(pr["moments"][:kmax])
        for s, v in acc.items():
            out.setdefault(s, np.full((len(rows), kmax), np.nan))[i] = np.mean(v, axis=0)
    return out


def liberation_curve(a, kmax=4, tmax=12.0, npts=400):
    """Free liberation (free Jacobi process, Demni-Hamdi-Hmidi 2012 eq. 2) of a projection P against itself:
    J_t = P Y_t P Y_t^* P, tau(P) = a  (lambda = 1, theta = a), m_n(0) = 1.  Returns t and m[n-1](t) = T_2n(t), n = 1..kmax
    (normalized by tau(P), i.e. our tau = Tr / d).  The Hilbert-space-isotropic null model (D -> infinity)."""
    from scipy.integrate import solve_ivp

    def rhs(t, m):
        mm = np.concatenate([[1.0], m])
        out = np.empty(kmax)
        for n in range(1, kmax + 1):
            acc = sum(mm[n - k - 1] * (mm[k] - mm[k + 1]) for k in range(0, n - 1))
            out[n - 1] = -n * mm[n] + a * n * mm[n - 1] + a * n * acc
        return out
    ts = np.linspace(0, tmax, npts)
    sol = solve_ivp(rhs, (0, tmax), np.ones(kmax), t_eval=ts, rtol=1e-10, atol=1e-12)
    return ts, sol.y


def chaos_fit(y, Tm, se, deg):
    """Weighted LSQ of T(y) = sum_{n<=deg} c_n y^n with T(1) = 1 imposed (eliminate c_deg)."""
    A = np.vander(y, deg + 1, increasing=True)
    w = 1.0 / np.maximum(se, 1e-12)
    # constraint sum c_n = 1  ->  c_deg = 1 - sum_{n<deg} c_n
    B = A[:, :-1] - A[:, -1:]
    b = Tm - A[:, -1]
    c, *_ = np.linalg.lstsq(B * w[:, None], b * w, rcond=None)
    c = np.append(c, 1.0 - c.sum())
    chi2 = float((((A @ c - Tm) * w) ** 2).sum())
    return c, chi2


def matched(T2_curve, Tk_curve, xs, T2_dec):
    """Interpolate the OU curve (parametrized by T2, monotone) at T2_dec; also return x_* (interpolated in x)."""
    order = np.argsort(T2_curve)
    t2, tk, x = T2_curve[order], Tk_curve[order], xs[order]
    if not np.all(np.diff(t2) > 0):
        return np.nan, np.nan
    return float(PchipInterpolator(t2, tk)(T2_dec)), float(PchipInterpolator(t2, x)(T2_dec))


def analyze(N, kmax=4, deg=5, nboot=400, seed=1):
    rows = load(N)
    a = rows[0]["a"]
    p = rows[0]["p"]
    xs, T = ou_table(rows, kmax)
    pr = probe_table(rows, kmax)
    ns = len(rows)
    rng = np.random.default_rng(seed)

    def stats(idx):
        Tm = T[idx].mean(0)
        c0 = Tm[-1, 1]                                  # T4 at x = 0 (last lag class, L = K/2)
        res = dict(T2=Tm[:, 0], T4=Tm[:, 1], R=Tm[:, 1] - c0 - 2 * a * (Tm[:, 0] - a), c0=c0)
        dec = {}
        for s, arr in pr.items():
            dm = np.nanmean(arr[idx], 0)
            m = {}
            for k in range(1, kmax):
                m[f"T{2 * k + 2}_ou"], m["x_star"] = matched(Tm[:, 0], Tm[:, k], xs, dm[0])
                m[f"T{2 * k + 2}_dec"] = dm[k]
            m["T2_dec"] = dm[0]
            m["T4_lib"] = matched(lib[0][::-1], lib[1][::-1], lib_t[::-1], dm[0])[0]
            m["R_dec"] = dm[1] - c0 - 2 * a * (dm[0] - a)
            m["R_ou"] = m["T4_ou"] - c0 - 2 * a * (dm[0] - a)
            dec[s] = m
        res["dec"] = dec
        return res

    lib_t, lib = liberation_curve(a, kmax)
    full = stats(np.arange(ns))
    boots = [stats(rng.integers(0, ns, ns)) for _ in range(nboot)]

    def se(f):
        v = np.array([f(b) for b in boots], dtype=float)
        return np.nanstd(v, axis=0, ddof=1)

    y = xs ** 2
    seT2, seT4 = se(lambda b: b["T2"]), se(lambda b: b["T4"])
    seT2[0] = seT4[0] = 1e-9
    fits = {}
    for name, Tm, s_ in (("T2", full["T2"], seT2), ("T4", full["T4"], seT4)):
        c, chi2 = chaos_fit(y, Tm, s_, deg)
        cb = np.array([chaos_fit(y, b[name], s_, deg)[0] for b in boots])
        fits[name] = dict(coef=c.tolist(), se=cb.std(0, ddof=1).tolist(), chi2=chi2, dof=len(y) - deg)
    out = dict(N=N, p=p, a=a, samples=ns, K=rows[0]["K"], x=xs.tolist(),
               T2=full["T2"].tolist(), T2_se=seT2.tolist(), T4=full["T4"].tolist(), T4_se=seT4.tolist(),
               R_ou=full["R"].tolist(), R_ou_se=se(lambda b: b["R"]).tolist(), c0=full["c0"],
               c0_se=float(se(lambda b: b["c0"])), chaos=fits, decoder={})
    for s, m in sorted(full["dec"].items()):
        e = {k: float(v) for k, v in m.items()}
        for k in m:
            e[k + "_se"] = float(se(lambda b, s=s, k=k: b["dec"][s][k]))
        e["x_exact"] = krawtchouk_x(N, p, s)
        e["x_ds"] = float(np.exp(-2 * p * s / N))
        e["dT4"] = e["T4_dec"] - e["T4_ou"]
        # three-way split of the excess over the fully free (Wachter) value 2a^2 - a^3
        e["T4_free"] = 2 * a * a - a ** 3
        e["split_lib"] = e["T4_lib"] - e["T4_free"]          # free liberation at the decoder's T_2 (analytic)
        e["split_iso"] = e["T4_ou"] - e["T4_lib"]            # isotropic SYK non-freeness (OU - liberation)
        e["split_loc"] = e["T4_dec"] - e["T4_ou"]            # localization (decoder - OU at equal T_2)
        e["dT4_se"] = float(se(lambda b, s=s: b["dec"][s]["T4_dec"] - b["dec"][s]["T4_ou"]))
        out["decoder"][int(s)] = e
    return out


def site_vs_ou(N, kmax=4, nboot=300, seed=2):
    """Site-local family T(1, z) (only couplings containing one site rotate) vs the isotropic OU curve, compared at
    equal T_2: excess_k(z) = T_2k^site(z) - T_2k^OU(at the same T_2).  z = -1 is the decoder of that site."""
    ou, st = load(N), load(N, site=N - 1)
    a = ou[0]["a"]
    xs, T = ou_table(ou, kmax)
    zs, Ts = ou_table(st, kmax)
    rng = np.random.default_rng(seed)

    def ex(io, js):
        To, Tz = T[io].mean(0), Ts[js].mean(0)
        c0 = To[-1, 1]
        e = {}
        for k in range(1, kmax):
            e[k] = np.array([matched(To[:, 0], To[:, k], xs, t2)[0] for t2 in Tz[:, 0]])
        xstar = np.array([matched(To[:, 0], To[:, 1], xs, t2)[1] for t2 in Tz[:, 0]])
        R = Tz[:, 1] - c0 - 2 * a * (Tz[:, 0] - a)
        return Tz, e, xstar, R, c0
    Tz, e, xstar, R, c0 = ex(np.arange(len(ou)), np.arange(len(st)))
    bs = [ex(rng.integers(0, len(ou), len(ou)), rng.integers(0, len(st), len(st))) for _ in range(nboot)]
    se4 = np.nanstd([b[0][:, 1] - b[1][1] for b in bs], axis=0, ddof=1)
    print(f"\n=== site family vs OU, N'={N}  a={a:.4f}  OU samples={len(ou)}  site samples={len(st)}")
    print("  z       T2_site  x*(OU)  T4_site  T4_OU(T2)  excess4           R_site   T6 excess  T8 excess")
    for i, z in enumerate(zs):
        if i == 0:
            continue
        print(f"  {z:+.3f}  {Tz[i, 0]:.5f}  {xstar[i]:.4f}  {Tz[i, 1]:.5f}  {e[1][i]:.5f}    "
              f"{Tz[i, 1] - e[1][i]:+.5f}±{se4[i]:.5f}  {R[i]:+.5f}  {Tz[i, 2] - e[2][i]:+.5f}   {Tz[i, 3] - e[3][i]:+.5f}")
    return dict(N=N, z=zs.tolist(), T2=Tz[:, 0].tolist(), T4=Tz[:, 1].tolist(), T4_ou_matched=e[1].tolist(),
                excess4_se=se4.tolist(), x_star=xstar.tolist(), R=R.tolist())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Np", type=int, nargs="+", required=True)
    ap.add_argument("--deg", type=int, default=5)
    ap.add_argument("--nboot", type=int, default=400)
    ap.add_argument("--json", default=None)
    ap.add_argument("--site", action="store_true", help="also compare the site-local family with the OU curve")
    ap.add_argument("--markdown", action="store_true", help="print the decoder summary table in markdown")
    a = ap.parse_args()
    allres = []
    for N in a.Np:
        try:
            r = analyze(N, deg=a.deg, nboot=a.nboot)
        except FileNotFoundError:
            print(f"N'={N}: no data")
            continue
        allres.append(r)
        print(f"\n=== N'={N}  a={r['a']:.4f}  samples={r['samples']}  c0=T4(0)={r['c0']:.5f}±{r['c0_se']:.5f}")
        print("  x       T2               T4               R_OU")
        for x, t2, s2, t4, s4, R, sR in zip(r["x"], r["T2"], r["T2_se"], r["T4"], r["T4_se"], r["R_ou"], r["R_ou_se"]):
            print(f"  {x:.3f}  {t2:.5f}±{s2:.5f}  {t4:.5f}±{s4:.5f}  {R:+.5f}±{sR:.5f}")
        for name in ("T2", "T4"):
            f = r["chaos"][name]
            print(f"  chaos {name}: " + " ".join(f"{c:+.4f}({e:.4f})" for c, e in zip(f["coef"], f["se"]))
                  + f"  chi2={f['chi2']:.1f}/{f['dof']}")
        print("  s  x_exact x_ds  x_*            T2_dec   T4_dec   T4_OU(x*)  dT4              R_dec    R_OU(x*)  "
              "T6_dec  T6_OU   T8_dec  T8_OU   T4_lib")
        for s, e in r["decoder"].items():
            print(f"  {s:<2d} {e['x_exact']:+.3f}  {e['x_ds']:.3f} {e['x_star']:.4f}±{e['x_star_se']:.4f} "
                  f"{e['T2_dec']:.5f}  {e['T4_dec']:.5f}  {e['T4_ou']:.5f}   {e['dT4']:+.5f}±{e['dT4_se']:.5f}  "
                  f"{e['R_dec']:+.5f}  {e['R_ou']:+.5f}   {e['T6_dec']:.5f} {e['T6_ou']:.5f}  "
                  f"{e['T8_dec']:.5f} {e['T8_ou']:.5f}  {e['T4_lib']:.5f}")
    print("\nSplit of T4_dec - T4_free (s = 1):  liberation | isotropic (OU - lib) | localization (dec - OU)")
    for r in allres:
        if "decoder" in r and 1 in r["decoder"]:
            e = r["decoder"][1]
            tot = e["T4_dec"] - e["T4_free"]
            print(f"  N'={r['N']:<2d} total {tot:.4f}:  {e['split_lib']:.4f} ({100 * e['split_lib'] / tot:.0f}%) | "
                  f"{e['split_iso']:+.4f} ({100 * e['split_iso'] / tot:.1f}%) | {e['split_loc']:.4f} "
                  f"({100 * e['split_loc'] / tot:.0f}%)")
    if a.site:
        for N in a.Np:
            try:
                allres.append(dict(site_family=site_vs_ou(N)))
            except FileNotFoundError:
                print(f"N'={N}: no site-family data")
    if a.markdown:
        print("\n| N′ | a | samples | T₂^dec | x_* | x_exact | T₄^dec | T₄^OU(x_*) | T₄^lib | T₄^free | liberation | isotropic | localization | R_dec | R_OU(x_*) |")
        print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for r in allres:
            if "decoder" not in r or 1 not in r["decoder"]:
                continue
            e = r["decoder"][1]
            tot = e["T4_dec"] - e["T4_free"]
            print(f"| {r['N']} | {r['a']:.4f} | {r['samples']} | {e['T2_dec']:.4f}({1e4 * e['T2_dec_se']:.0f}) | "
                  f"{e['x_star']:.3f} | {e['x_exact']:.3f} | {e['T4_dec']:.4f}({1e4 * e['T4_dec_se']:.0f}) | "
                  f"{e['T4_ou']:.4f}({1e4 * e['T4_ou_se']:.0f}) | {e['T4_lib']:.4f} | {e['T4_free']:.4f} | "
                  f"{e['split_lib']:.4f} ({100 * e['split_lib'] / tot:.0f} %) | {e['split_iso']:+.4f} "
                  f"({100 * e['split_iso'] / tot:.1f} %) | {e['split_loc']:.4f}({1e4 * e['dT4_se']:.0f}) "
                  f"({100 * e['split_loc'] / tot:.0f} %) | {e['R_dec']:+.4f} | {e['R_ou']:+.4f} |")
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(allres, fh, indent=1)


if __name__ == "__main__":
    main()
