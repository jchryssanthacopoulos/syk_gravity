#!/usr/bin/env python3
r"""
Reproduce Figure 12: adjacent-gap statistics of the nontrivial down-channel
decoder spectrum in the one-flavor N=2 SYK model (the generic q=3 wedge model).

Setup (Section 8.4).  The down-channel decoder operator acts entirely inside the
enlarged BPS space B_{N'}^P:

    O_down = P_B (1 - n_{N'}) P_B ,

a mode-occupation projector (added mode empty) compressed to the highly
degenerate BPS subspace -- an LMRS-type operator.  Its eigenvalues lambda in [0,1]
are the squared decoder singular values.  We strip the structural endpoints
lambda = 0 (states with the added mode fully occupied) and lambda = 1 (added mode
fully empty), and study the level statistics of the interior spectrum via the
adjacent-gap ratio (Eq. 98)

    s_i = lambda_{i+1} - lambda_i ,   r_i = min(s_i, s_{i+1}) / max(s_i, s_{i+1}),

which needs no unfolding.  Pooled over independent complex-Gaussian realizations,
<r> is compared with Poisson / GOE / GUE (Eq. 99), and P(r), F(r) with Fig. 12.

Model machinery (gen_form, harmonic, subsets) is imported from reproduce_table1.py:
H = Lambda^*(C^{N'}), Q = C ^ (.) with a generic 3-form C, BPS = harmonic.  Only
the enlarged theory is needed; O_down does not reference the size-N couplings.

The paper's figure is N'=12, P=6, 100 realizations (Eq. 100: <r>=0.59853+/-0.00203).

    python decoder_chaos.py --nprime 12 --P 6 --realizations 100 --plot fig12.pdf
    python decoder_chaos.py --nprime 8  --P 4 --realizations 50           # quick test
"""
import argparse
import numpy as np
from reproduce_table1 import gen_form, harmonic, subsets

_trapz = getattr(np, "trapezoid", np.trapezoid)   # numpy>=2 renamed trapz

# ------------------------------------------------------------------ spectrum ---
def decoder_eigs(Nprime, P, seed):
    """Eigenvalues of O_down = P_B (1 - n_{N'}) P_B on the BPS space B_{N'}^P."""
    f = gen_form(Nprime, seed)
    B = harmonic(f, Nprime, P)                     # (C(N',P) x dimB), orthonormal cols
    if B.shape[1] == 0:
        return np.array([])
    subs = subsets(Nprime, P)
    empty = np.array([0.0 if (Nprime in S) else 1.0 for S in subs])   # (1 - n_{N'})
    O = B.conj().T @ (empty[:, None] * B)          # B^dagger diag(1-n) B
    O = (O + O.conj().T) / 2
    return np.linalg.eigvalsh(O).real

def interior_and_ratios(eigs, endtol=1e-6):
    """Strip structural lambda=0,1; return (interior eigenvalues, gap ratios r)."""
    interior = np.sort(eigs[(eigs > endtol) & (eigs < 1.0 - endtol)])
    if len(interior) < 3:
        return interior, np.array([])
    s = np.diff(interior)
    r = np.minimum(s[:-1], s[1:]) / np.maximum(s[:-1], s[1:])
    return interior, r

# --------------------------------------------------------- RMT surmise (r~) ----
_RMT_MEAN = {"Poisson": 0.38629, "GOE": 0.53590, "GUE": 0.60266}

def _surmise_unnorm(r, beta):
    if beta == 0:                                  # Poisson
        return 2.0 / (1.0 + r) ** 2
    return (r + r * r) ** beta / (1.0 + r + r * r) ** (1.0 + 1.5 * beta)

def rmt_pdf(r, kind):
    beta = {"Poisson": 0, "GOE": 1, "GUE": 2}[kind]
    grid = np.linspace(0, 1, 20001)
    Z = _trapz(_surmise_unnorm(grid, beta), grid)   # normalize on [0,1]
    return _surmise_unnorm(r, beta) / Z

def rmt_cdf(r, kind):
    beta = {"Poisson": 0, "GOE": 1, "GUE": 2}[kind]
    grid = np.linspace(0, 1, 20001)
    dens = _surmise_unnorm(grid, beta)
    dens /= _trapz(dens, grid)
    cdf = np.concatenate([[0], np.cumsum((dens[1:] + dens[:-1]) / 2 * np.diff(grid))])
    return np.interp(r, grid, cdf)

# ------------------------------------------------------------------- driver ----
def run(Nprime, P, realizations, seed0, endtol):
    all_r, per_real_mean, n_interior = [], [], []
    for k in range(realizations):
        eigs = decoder_eigs(Nprime, P, seed0 + k)
        interior, r = interior_and_ratios(eigs, endtol)
        n_interior.append(len(interior))
        if len(r):
            all_r.append(r)
            per_real_mean.append(r.mean())
        if (k + 1) % 25 == 0 or k == realizations - 1:
            print(f"  realization {k+1}/{realizations} done", flush=True)
    r = np.concatenate(all_r) if all_r else np.array([])
    per_real_mean = np.array(per_real_mean)
    return r, per_real_mean, np.array(n_interior)

def make_plot(r, Nprime, P, outfile):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    grid = np.linspace(1e-4, 1, 400)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 4.2))
    # left: density
    axL.hist(r, bins=40, range=(0, 1), density=True, color="0.8",
             edgecolor="0.5", label="data")
    styles = {"Poisson": "C7--", "GOE": "C0-", "GUE": "C3-"}
    for kind, st in styles.items():
        axL.plot(grid, rmt_pdf(grid, kind), st, lw=1.8, label=kind)
    axL.set_xlabel("$r$"); axL.set_ylabel("$P(r)$"); axL.set_xlim(0, 1)
    axL.legend(); axL.set_title(f"decoder gap ratios  $N'={Nprime}$, $P={P}$")
    # right: CDF + KS
    rs = np.sort(r); Fe = np.arange(1, len(rs) + 1) / len(rs)
    axR.plot(rs, Fe, "k-", lw=1.6, label="data")
    ks = {}
    for kind, st in styles.items():
        Ft = rmt_cdf(rs, kind)
        ks[kind] = np.max(np.abs(Fe - Ft))
        axR.plot(grid, rmt_cdf(grid, kind), st, lw=1.5,
                 label=f"{kind} (KS={ks[kind]:.3f})")
    axR.set_xlabel("$r$"); axR.set_ylabel("$F(r)$"); axR.set_xlim(0, 1); axR.set_ylim(0, 1)
    axR.legend(loc="lower right"); axR.set_title("cumulative")
    fig.tight_layout(); fig.savefig(outfile, dpi=150)
    print(f"  wrote {outfile}   KS: " + ", ".join(f"{k}={v:.4f}" for k, v in ks.items()))
    return ks

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nprime", type=int, default=12, help="enlarged size N'")
    ap.add_argument("--P", type=int, default=6, help="total charge P (central: N'/2)")
    ap.add_argument("--realizations", type=int, default=100)
    ap.add_argument("--seed0", type=int, default=0)
    ap.add_argument("--endtol", type=float, default=1e-6,
                    help="tolerance for stripping structural lambda=0,1")
    ap.add_argument("--plot", default=None, help="output figure file (pdf/png)")
    a = ap.parse_args()

    print(f"one-flavor N=2 SYK decoder chaos:  N'={a.nprime}, P={a.P}, "
          f"{a.realizations} realizations")
    r, per_real_mean, n_int = run(a.nprime, a.P, a.realizations, a.seed0, a.endtol)
    if len(r) == 0:
        print("no interior eigenvalues -- try a larger/central sector"); return
    mean_r = per_real_mean.mean()
    se = per_real_mean.std(ddof=1) / np.sqrt(len(per_real_mean)) if len(per_real_mean) > 1 else 0.0
    print(f"\n  interior eigenvalues / realization: {n_int.mean():.1f} "
          f"(total {n_int.sum()}),  gap ratios: {len(r)}")
    print(f"  <r> = {mean_r:.5f} +/- {se:.5f}")
    print("  benchmarks:  " + "  ".join(f"{k} {v}" for k, v in _RMT_MEAN.items()))
    closest = min(_RMT_MEAN, key=lambda k: abs(_RMT_MEAN[k] - mean_r))
    print(f"  -> closest benchmark: {closest}")
    if a.plot:
        make_plot(r, a.nprime, a.P, a.plot)

if __name__ == "__main__":
    main()
