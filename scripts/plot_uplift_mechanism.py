#!/usr/bin/env python3
r"""
Mechanism figure for the heterogeneity prediction (docs/derivations.md D6 §§5-6).

Everything is computed from the OLD model (N'-1 modes) through the exact uplift (src/uplift.py), in the space
X (+) Y = Lambda^F(o) (+) Lambda^{F-1}(o) ~ Lambda^F(N'):
  first order:   y1 = ((1 - s^2)/(1 + s^2))^2 for s = sigma (h-states) and s = sigma' (g-states);
  self-dressed:  spectra of the diagonal blocks M_h, M_g (h-g coupling dropped);
  exact:         full block operator (= exact N'-mode decoder, tests/test_uplift.py).
Per-state tracking: first-order states Psi1 = (h_j, -eps K h_j), (eps K' g_j, g_j) (singular vectors), projected onto
the exact BPS space and Loewdin-orthonormalized; T[i, j] = |<exact_i | psi_j>|^2 is doubly stochastic, so
mean_j E[y | j] = T_2(exact) exactly.

Panels: (a) distribution of y (first / self-dressed / exact); (b) per-state shifts E[y|j] - y1 and
E[y(1-y)|j] - y1(1-y1), binned in y1.

Usage: .venv/bin/python scripts/plot_uplift_mechanism.py --Np 14 --seeds 9000 9001
Memory: N' = 14 ~ 0.5 GB; N' = 15 ~ 2 GB.
"""
import argparse
import json
import os
import sys

import numpy as np
import scipy.linalg as sla

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from chord_chaos import Sector, random_couplings  # noqa: E402
from uplift import OldModel  # noqa: E402
from run_uplift_walls import git_state  # noqa: E402


def one_seed(Np, seed):
    rng = np.random.default_rng(seed)
    c = random_couplings(rng, Sector(Np, 3, F=Np // 2).M)
    om = OldModel(Np, c, rng)
    eps = om.eps
    a, b = om.uplift_components()
    W = np.vstack([a, b])                                           # exact BPS states (columns), X (+) Y
    G = W.conj().T @ W
    N = a.conj().T @ a - b.conj().T @ b
    mu, V = sla.eigh(N, G)                                          # V^dag G V = 1 -> W V orthonormal
    X = W @ V
    nh = om.VF.shape[1]
    mu_self = np.concatenate([sla.eigvalsh(N[:nh, :nh], G[:nh, :nh]), sla.eigvalsh(N[nh:, nh:], G[nh:, nh:])])
    KH, KpH, _ = om.three_operators()
    U1, s1, W1h = np.linalg.svd(KH, full_matrices=True)
    sig = np.zeros(nh); sig[:len(s1)] = s1
    h = om.VF @ W1h.conj().T
    U2, s2, W2h = np.linalg.svd(KpH, full_matrices=True)
    sigp = np.zeros(om.VF1.shape[1]); sigp[:len(s2)] = s2
    g = om.VF1 @ W2h.conj().T
    Psi1 = np.vstack([np.hstack([h, eps * (KpH @ W2h.conj().T)]), np.hstack([-eps * (KH @ W1h.conj().T), g])])
    Psi1 /= np.linalg.norm(Psi1, axis=0)
    mu1 = np.concatenate([(1 - sig ** 2) / (1 + sig ** 2), -(1 - sigp ** 2) / (1 + sigp ** 2)])
    C = X.conj().T @ Psi1
    Uc, _, Vch = np.linalg.svd(C)
    T = np.abs(Uc @ Vch) ** 2                                       # doubly stochastic transfer matrix
    y, y1 = mu ** 2, mu1 ** 2
    return dict(y=y, y_self=mu_self ** 2, y1=y1, Ey=T.T @ y, Ef=T.T @ (y * (1 - y)),
                capture=np.linalg.norm(C, axis=0) ** 2)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Np", type=int, default=14)
    ap.add_argument("--seeds", type=int, nargs="+", default=[9000, 9001])
    ap.add_argument("--out", default=os.path.join(ROOT, "results/figures/uplift_2026-10-12/fig_uplift_mechanism.pdf"))
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    runs = [one_seed(a.Np, s) for s in a.seeds]
    cat = {k: np.concatenate([r[k] for r in runs]) for k in runs[0]}
    H = lambda y: float(np.mean(y * (1 - y)))  # noqa: E731
    summ = dict(Np=a.Np, seeds=a.seeds, git=git_state(),
                T2=dict(first=float(cat["y1"].mean()), self=float(cat["y_self"].mean()), exact=float(cat["y"].mean())),
                H=dict(first=H(cat["y1"]), self=H(cat["y_self"]), exact=H(cat["y"])),
                capture_mean=float(cat["capture"].mean()), capture_min=float(cat["capture"].min()))
    bins = np.linspace(0, 1, 11)
    rows = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        s = (cat["y1"] >= lo) & ((cat["y1"] < hi) if hi < 1 else (cat["y1"] <= 1))
        if s.sum():
            rows.append(dict(lo=float(lo), hi=float(hi), n=int(s.sum()), dy=float(np.mean(cat["Ey"][s] - cat["y1"][s])),
                             dH=float(np.mean(cat["Ef"][s] - cat["y1"][s] * (1 - cat["y1"][s])))))
    summ["bins"] = rows
    print(json.dumps({k: v for k, v in summ.items() if k != "bins"}, indent=1))
    for r in rows:
        print(f"  y1 in [{r['lo']:.1f},{r['hi']:.1f}): n={r['n']:5d}  E[y|j]-y1 = {r['dy']:+.4f}   "
              f"E[y(1-y)|j]-y1(1-y1) = {r['dH']:+.4f}")
    if a.json:
        json.dump(summ, open(a.json, "w"), indent=1)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.labelsize": 10, "legend.fontsize": 8.5})
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.6))
    hb = np.linspace(0, 1, 26)
    for key, lab, sty in (("y1", "first order", dict(color="tab:red", ls="--")),
                          ("y_self", "self-dressed (blocks)", dict(color="tab:blue", ls=":")),
                          ("y", "exact", dict(color="k", ls="-"))):
        yy = cat[key]
        ax[0].hist(yy, bins=hb, density=True, histtype="step", lw=1.4,
                   label=f"{lab}: $T_2$={yy.mean():.3f}, $H$={H(yy):.4f}", **sty)
    ax[0].set_yscale("log")
    from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
    ax[0].yaxis.set_major_locator(FixedLocator([0.3, 0.5, 1, 2, 5]))
    ax[0].yaxis.set_major_formatter(FixedFormatter(["0.3", "0.5", "1", "2", "5"]))
    ax[0].yaxis.set_minor_locator(NullLocator())
    ax[0].set_xlabel(r"$y=\mu^2=\cos^2 2\theta$")
    ax[0].set_ylabel("density")
    ax[0].set_title(rf"decoder angles, $N'={a.Np}$ ({len(a.seeds)} realizations)", fontsize=10)
    ax[0].legend(frameon=False, loc="upper center")
    xc = np.array([(r["lo"] + r["hi"]) / 2 for r in rows])
    ax[1].plot(xc, [r["dy"] for r in rows], "o-", color="tab:purple", label=r"$E[y\,|\,j]-y_1$")
    ax[1].plot(xc, [r["dH"] for r in rows], "s-", color="tab:green", mfc="none",
               label=r"$E[y(1-y)\,|\,j]-y_1(1-y_1)$")
    ax[1].axhline(0, color="k", lw=0.6)
    ax[1].set_xlabel(r"first-order $y_1$ of state $j$")
    ax[1].set_ylabel("per-state change (exact $-$ first order)")
    ax[1].set_title("first order $\\to$ exact: shift, not reflection", fontsize=10)
    ax[1].legend(frameon=False, loc="lower left")
    fig.tight_layout()
    fig.savefig(a.out)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
