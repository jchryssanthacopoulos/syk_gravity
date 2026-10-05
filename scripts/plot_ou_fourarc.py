#!/usr/bin/env python3
r"""
Figure for Phase A1/A2 (docs/derivations.md D5): the decoder fourth moment against the isotropic OU family.

(a) The (T_2, T_4) plane: OU curves (isotropic coupling rotation, T(x, x)), the site-local family T(1, z) (only the
    couplings containing one site rotate; z = -1 is the decoder), the decoder (s = 1) and parity probes (s >= 2).
(b) The remainder R = T_4 - c0 - 2a(T_2 - a) (c0 = two-model value) against eta = (T_2 - a)/(1 - a).
(c) Excess of T_4 over the OU curve at equal T_2, vs N': decoder (s = 1) and the site family at z = 0 (Gaussian
    resampling of the site's couplings).
(d) Principal-angle spectra (cos^2 theta) at N'_max: decoder vs the OU pair whose T_2 is matched (histograms linearly
    interpolated in T_2 between the two neighbouring lags).

Usage: .venv/bin/python scripts/plot_ou_fourarc.py --Np 8 9 10 11 12 13 14 --out results/figures/ou_fourarc_2026-10-05
"""
import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_ou_fourarc import liberation_curve, load, matched, ou_table, probe_table  # noqa: E402


def liberation_spectrum(a, t, D=2000, steps=40, seed=11):
    """Finite-D free liberation: P = diag(1_{aD}, 0), Y_t = prod exp(i sqrt(dt) H_k) with GUE H_k (tau(H^2) = 1), so
    tau(Y_t) -> exp(-t/2) (free unitary Brownian motion).  cos^2 of the principal angles between P and Y P Y^* are the
    squared singular values of the aD x aD block of Y_t."""
    import scipy.linalg as sla
    rng = np.random.default_rng(seed)
    k = int(round(a * D))
    Y = np.eye(D, dtype=complex)
    dt = t / steps
    for _ in range(steps):
        X = (rng.standard_normal((D, D)) + 1j * rng.standard_normal((D, D))) / np.sqrt(2)
        H = (X + X.conj().T) / np.sqrt(2 * D)
        w, U = sla.eigh(H)
        Y = (U * np.exp(1j * np.sqrt(dt) * w)) @ (U.conj().T @ Y)
    return sla.svdvals(Y[:k, :k]) ** 2


def curve(rows, kmax=4):
    xs, T = ou_table(rows, kmax)
    return xs, T.mean(0), T.std(0, ddof=1) / np.sqrt(len(rows))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Np", type=int, nargs="+", required=True)
    ap.add_argument("--panel-a", type=int, nargs="+", default=[12, 14])
    ap.add_argument("--out", required=True)
    ap.add_argument("--lib-sim", type=int, default=0, help="finite-D liberation spectrum in panel (d) (D, e.g. 2000)")
    ap.add_argument("--figsize", type=float, nargs=2, default=[9.2, 7.6])
    ap.add_argument("--fontsize", type=float, default=9)
    ap.add_argument("--name", default="ou_fourarc", help="output basename")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    fs = a.fontsize
    plt.rcParams.update({"font.size": fs, "axes.labelsize": fs + 1, "legend.fontsize": fs - 1.5,
                         "axes.titlesize": fs, "xtick.labelsize": fs - 0.5, "ytick.labelsize": fs - 0.5})
    fig, ax = plt.subplots(2, 2, figsize=tuple(a.figsize))
    cols = plt.cm.viridis(np.linspace(0.0, 0.85, len(a.Np)))
    summary = []
    for N, col in zip(a.Np, cols):
        try:
            ou = load(N)
        except FileNotFoundError:
            continue
        aa = ou[0]["a"]
        xs, Tm, Te = curve(ou)
        c0 = Tm[-1, 1]
        pr = probe_table(ou, 4)
        dec = {s: np.nanmean(v, 0) for s, v in pr.items()}
        dec_se = {s: np.nanstd(v, 0, ddof=1) / np.sqrt(np.isfinite(v[:, 0]).sum()) for s, v in pr.items()}
        try:
            st = load(N, site=N - 1)
            zs, Tz, _ = curve(st)
        except FileNotFoundError:
            st = None
        # (a)
        if N in a.panel_a:
            ax[0, 0].plot(Tm[:, 0], Tm[:, 1], "-o", color=col, ms=3, lw=1.2, label=f"OU $T(x,x)$, $N'={N}$")
            if st is not None:
                ax[0, 0].plot(Tz[:, 0], Tz[:, 1], "--s", color=col, ms=3, lw=1, mfc="none",
                              label=f"site-local $T(1,z)$, $N'={N}$")
            for s, m in sorted(dec.items()):
                ax[0, 0].errorbar(m[0], m[1], xerr=dec_se[s][0], yerr=dec_se[s][1], fmt="*" if s == 1 else "^",
                                  color=col, ms=9 if s == 1 else 5, mec="k", mew=0.5,
                                  label=(f"decoder $s=1$ / parity $s\\geq2$, $N'={N}$" if s == 1 else None))
        # (b)
        eta = (Tm[:, 0] - aa) / (1 - aa)
        R = Tm[:, 1] - c0 - 2 * aa * (Tm[:, 0] - aa)
        ax[0, 1].plot(eta, R, "-", color=col, lw=1.3, label=f"OU, $N'={N}$")
        if st is not None:
            etaz = (Tz[:, 0] - aa) / (1 - aa)
            Rz = Tz[:, 1] - c0 - 2 * aa * (Tz[:, 0] - aa)
            ax[0, 1].plot(etaz, Rz, "--", color=col, lw=1)
        for s, m in sorted(dec.items()):
            ax[0, 1].plot((m[0] - aa) / (1 - aa), m[1] - c0 - 2 * aa * (m[0] - aa), "*" if s == 1 else "^",
                          color=col, ms=9 if s == 1 else 4, mec="k", mew=0.5)
        # (c)
        e1 = dec[1][1] - matched(Tm[:, 0], Tm[:, 1], xs, dec[1][0])[0]
        ez = None
        if st is not None:
            iz = int(np.argmin(np.abs(zs)))
            ez = Tz[iz, 1] - matched(Tm[:, 0], Tm[:, 1], xs, Tz[iz, 0])[0]
        lt, lm = liberation_curve(aa, 2)
        t4lib = matched(lm[0][::-1], lm[1][::-1], lt[::-1], dec[1][0])[0]
        t4ou = matched(Tm[:, 0], Tm[:, 1], xs, dec[1][0])[0]
        free = 2 * aa * aa - aa ** 3
        if N in a.panel_a:
            ax[0, 0].plot(lm[0], lm[1], ":", color=col, lw=1, label=f"free liberation, $a={aa:.3f}$")
        summary.append(dict(N=N, a=aa, dec_excess=e1, site_z0_excess=ez, split_lib=t4lib - free,
                            split_iso=t4ou - t4lib, split_loc=dec[1][1] - t4ou, T4_dec=dec[1][1], T4_free=free,
                            R1=1 - c0 - 2 * aa * (1 - aa), R_dec=dec[1][1] - c0 - 2 * aa * (dec[1][0] - aa),
                            R_ou_matched=matched(Tm[:, 0], Tm[:, 1], xs, dec[1][0])[0] - c0 - 2 * aa * (dec[1][0] - aa)))
    # (c) stacked split of T4_dec - T4_free
    Ns = [r["N"] for r in summary]
    lib = np.array([r["split_lib"] for r in summary])
    iso = np.array([r["split_iso"] for r in summary])
    loc = np.array([r["split_loc"] for r in summary])
    tot = lib + iso + loc
    ax[1, 0].bar(Ns, lib / tot, color="0.75", label=r"free liberation at $T_2^{\rm dec}$ (analytic)")
    ax[1, 0].bar(Ns, iso / tot, bottom=lib / tot, color="C0", label="isotropic non-freeness (OU $-$ lib.)")
    ax[1, 0].bar(Ns, loc / tot, bottom=(lib + iso) / tot, color="C3", label="localization (decoder $-$ OU)")
    zz = [(r["N"], (r["split_lib"] + r["split_iso"] + r["site_z0_excess"]) / (r["T4_dec"] - r["T4_free"]))
          for r in summary if r["site_z0_excess"] is not None]
    if zz:
        ax[1, 0].plot(*zip(*zz), "k_", ms=14, mew=2, label="part of red reached by site resampling ($z=0$)")
    ax[1, 0].set_xlabel("$N'$")
    ax[1, 0].set_ylabel(r"share of $T_4^{\rm dec}-T_4^{\rm free}$")
    ax[1, 0].set_ylim(0, 1.3)
    ax[1, 0].legend(fontsize=fs - 2, loc="upper left", ncol=1)
    ax[1, 0].set_title("(c) decoder fourth-moment excess over Wachter, by origin", fontsize=fs)
    # (d) spectra at the largest N'
    Nmax = max(Ns)
    ou = load(Nmax)
    bins = len(ou[0]["ou"][0]["hist"])
    edges = np.linspace(0, 1, bins + 1)
    cen = 0.5 * (edges[1:] + edges[:-1])
    H = {o["lag"]: np.zeros(bins) for o in ou[0]["ou"]}
    T2l = {o["lag"]: [] for o in ou[0]["ou"]}
    for r in ou:
        for o in r["ou"]:
            H[o["lag"]] += np.array(o["hist"], float)
            T2l[o["lag"]].append(o["moments"][0])
    hd = np.zeros(bins)
    t2d = []
    for r in ou:
        for p_ in r["probes"]:
            if p_["s"] == 1:
                hd += np.array(p_["hist"], float)
                t2d.append(p_["moments"][0])
    t2d = np.mean(t2d)
    lags = sorted(H)
    t2s = np.array([np.mean(T2l[L]) for L in lags])
    j = int(np.searchsorted(-t2s, -t2d))                 # t2s decreasing with lag
    L1, L2 = lags[j - 1], lags[j]
    w = (t2d - t2s[j]) / (t2s[j - 1] - t2s[j])
    hou = w * H[L1] / H[L1].sum() + (1 - w) * H[L2] / H[L2].sum()
    width = edges[1] - edges[0]
    ax[1, 1].step(cen, hd / hd.sum() / width, where="mid", color="C3", lw=1.2, label=f"decoder, $T_2={t2d:.3f}$")
    ax[1, 1].step(cen, hou / width, where="mid", color="C0", lw=1.2, label="OU at matched $T_2$ (interp.)")
    if a.lib_sim > 0:
        aN = ou[0]["a"]
        tl = -np.log((t2d - aN) / (1 - aN))
        c2l = liberation_spectrum(aN, tl, D=a.lib_sim)
        nb = 50                                           # coarser bins: ~aD eigenvalues only
        hl, el = np.histogram(c2l, bins=nb, range=(0, 1))
        ax[1, 1].step(0.5 * (el[1:] + el[:-1]), hl / hl.sum() / (1.0 / nb), where="mid", color="0.4", lw=1.0, ls="--",
                      label=f"free liberation, $t={tl:.2f}$ ($D={a.lib_sim}$, $T_2={c2l.mean():.3f}$)")
    ax[1, 1].set_yscale("log")
    ax[1, 1].set_xlabel(r"$\cos^2\theta$ (principal angles of $P$, $P'$)")
    ax[1, 1].set_ylabel("density")
    ax[1, 1].set_title(f"(d) principal-angle spectra, $N'={Nmax}$", fontsize=fs)
    ax[1, 1].legend()
    ax[0, 0].set_xlabel(r"$T_2=\tau(PP')$")
    ax[0, 0].set_ylabel(r"$T_4=\tau(PP'PP')$")
    ax[0, 0].set_title("(a) fourth vs second moment", fontsize=fs)
    ax[0, 0].legend(ncol=1)
    ax[0, 1].axhline(0, color="k", lw=0.5)
    ax[0, 1].set_xlabel(r"$\eta=(T_2-a)/(1-a)$")
    ax[0, 1].set_ylabel(r"$R=T_4-c_0-2a(T_2-a)$")
    ax[0, 1].set_title("(b) remainder $R$: OU (—), site-local (- -), decoder ($\\star$)", fontsize=fs)
    ax[0, 1].legend()
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(a.out, f"{a.name}.{ext}"), dpi=180)
    with open(os.path.join(a.out, "ou_fourarc_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    for r in summary:
        print(r)


if __name__ == "__main__":
    main()
