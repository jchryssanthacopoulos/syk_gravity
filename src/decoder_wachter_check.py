#!/usr/bin/env python3
r"""
Exact-diagonalization decoder spectrum vs the analytic Wachter prediction.

The Wachter parameters are taken from CLOSED-FORM dimensions (no fitting):
  ambient   D = C(N+1, p)
  slot      m = C(N, p)                      -> b = m/D = (N+1-p)/(N+1)
  BPS dim   d = Dtil_p - Dtil_{N+1-p-q}      -> a = d/D     (index/walk formula)
and the interior density is
  f(l) = sqrt((l_+ - l)(l - l_-)) / (2 pi a l (1-l)),   l_+- = (sqrt(a(1-b)) +- sqrt(b(1-a)))^2 .

Claim tested: the *interior* (0<l<1) exact eigenvalues follow this curve, while the
endpoints l=0,1 carry protected excess beyond the Wachter atoms.

    python decoder_wachter_check.py --nmin 6 --nmax 10 --realizations 60 --plot decoder_wachter.png
"""
import argparse
from math import comb, pi, sqrt
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from decoder_chaos import decoder_eigs
q = 3
TOL = 1e-6
_trapz = getattr(np, "trapezoid", np.trapezoid)


def Dtil(M, k):
    s, n = 0, 0
    while k - n * q >= 0:
        s += (-1) ** n * comb(M, k - n * q); n += 1
    return s


def wachter(a, b, grid):
    lm = (sqrt(a * (1 - b)) - sqrt(b * (1 - a))) ** 2
    lp = (sqrt(a * (1 - b)) + sqrt(b * (1 - a))) ** 2
    f = np.zeros_like(grid)
    inb = (grid > lm) & (grid < lp)
    f[inb] = np.sqrt(np.clip((lp - grid[inb]) * (grid[inb] - lm), 0, None)) \
        / (2 * pi * a * grid[inb] * (1 - grid[inb]))
    p0 = max(0.0, 1 - b / a); p1 = max(0.0, 1 - (1 - b) / a)
    return f, p0, p1, lm, lp


def wachter_cdf(a, b, x):                       # CDF of the interior (unit-mass) law
    grid = np.linspace(1e-6, 1 - 1e-6, 6000)
    f, p0, p1, lm, lp = wachter(a, b, grid)
    area = _trapz(f, grid)
    cdf = np.concatenate([[0], np.cumsum((f[1:] + f[:-1]) / 2 * np.diff(grid))]) / area
    return np.interp(x, grid, cdf)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nmin", type=int, default=6)
    ap.add_argument("--nmax", type=int, default=10)
    ap.add_argument("--realizations", type=int, default=60)
    ap.add_argument("--plot", default="decoder_wachter.png")
    a_ = ap.parse_args()
    Ns = list(range(a_.nmin, a_.nmax + 1))

    ncol = 2; nrow = int(np.ceil(len(Ns) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 3.7 * nrow), squeeze=False)
    bins = np.linspace(0, 1, 33); ctr = 0.5 * (bins[1:] + bins[:-1])
    grid = np.linspace(1e-4, 1 - 1e-4, 1000)

    for idx, N in enumerate(Ns):
        ax = axes[idx // ncol][idx % ncol]
        Np, P = N + 1, (N + 1) // 2                 # half filling of enlarged theory
        D = comb(Np, P); m = comb(N, P)
        d = Dtil(Np, P) - Dtil(Np, Np - P - q)      # BPS dim from the closed-form index
        a = d / D; b = m / D                        # analytic Wachter parameters

        eigs = np.concatenate([decoder_eigs(Np, P, s) for s in range(a_.realizations)])
        inte = eigs[(eigs > TOL) & (eigs < 1 - TOL)]
        e0 = np.mean([np.mean(decoder_eigs(Np, P, s) < TOL) for s in range(a_.realizations)])
        e1 = np.mean([np.mean(decoder_eigs(Np, P, s) > 1 - TOL) for s in range(a_.realizations)])

        f, p0, p1, lm, lp = wachter(a, b, grid)
        area = _trapz(f, grid)
        ax.hist(inte, bins=bins, density=True, color="#7fb2e0", alpha=0.6,
                edgecolor="#33608f", label="exact diag (interior)")
        ax.plot(grid, f / area, "k-", lw=1.8, label="Wachter (closed form)")
        for e in (lm, lp):
            ax.axvline(e, ls=":", color="0.55", lw=1.0)
        # KS distance between exact interior ECDF and analytic Wachter CDF
        xs = np.sort(inte); Fe = np.arange(1, len(xs) + 1) / len(xs)
        ks = np.max(np.abs(Fe - wachter_cdf(a, b, xs))) if len(xs) else float("nan")
        ax.set_title(fr"${N}\!\to\!{N+1}$   $a={a:.3f},\ b={b:.3f}$", fontsize=10)
        ax.set_xlabel(r"$\lambda=\sigma^2$"); ax.set_xlim(0, 1)
        ax.text(0.03, 0.97,
                f"atoms $\\lambda{{=}}0$: exact {e0:.2f} / Wachter {p0:.2f}\n"
                f"atoms $\\lambda{{=}}1$: exact {e1:.2f} / Wachter {p1:.2f}\n"
                f"bulk KS = {ks:.3f}",
                transform=ax.transAxes, va="top", ha="left", fontsize=7.5,
                bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))
        if idx % ncol == 0:
            ax.set_ylabel("interior density")
        print(f"N={N}->{N+1}: a={a:.3f} b={b:.3f} edges=({lm:.3f},{lp:.3f})  "
              f"bulk KS={ks:.3f}  atom0 exact/Wachter={e0:.2f}/{p0:.2f}  "
              f"atom1 exact/Wachter={e1:.2f}/{p1:.2f}")
    for j in range(len(Ns), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    axes[0][-1].legend(fontsize=8, loc="upper right")
    fig.suptitle("Decoder spectrum: exact diagonalization vs analytic Wachter law "
                 "(one-flavor, down channel, half filling)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(a_.plot, dpi=160, bbox_inches="tight")
    print(f"wrote {a_.plot}")


if __name__ == "__main__":
    main()
