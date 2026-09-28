"""
Reproduction of Figure 5 of C. V. Johnson,
"Fortuitous Chaos, BPS Black Holes, and Random Matrices", arXiv:2601.17122, Sec. VI.

The SAME Wishart data as Fig. 4 (M is N x (N+Gamma) complex Gaussian, H = M^dag M with
Gamma exact zeros, lambda = lambda_tilde/(4N)), now ZOOMED into the hard edge by the
energy rescaling

    E = 4(Gamma_tilde + 2) N^2 lambda,      Gamma_tilde = Gamma/N,

which sets hbar = mu = 1 and moves the continuum edge to E_0 = Gamma^2. The finite-N
histogram then resolves into the double-scaled edge:

  leading (eq. 37):  rho_0(E) = (1/2 pi) sqrt(E - E_0) / E,          E > E_0,
  exact  (eq. 26):   rho(E)   = (1/4)[ J_G^2(xi) + J_{G+1}^2(xi)
                                       - (2G/xi) J_G(xi) J_{G+1}(xi) ],  xi = sqrt(E), G = Gamma.

rho(E) (Bessel kernel) captures the oscillations = the individual near-edge energy levels;
rho_0(E) is their smooth envelope. The Gamma exact zeros are the BPS states at E = 0.

Normalization note: the curves are a DENSITY OF STATES (levels per unit E, per sample), so
the histogram is normalized by counts/(reps * binwidth) -- NOT density=True.

Usage:  python3 reproduce_johnson_fig5.py       Requires: numpy, scipy, matplotlib.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv

# ---- parameters ----
N     = 100
Gamma = 10
reps  = 30000          # Johnson uses 1e5; 3e4 already resolves the first few level bumps
Emax  = 7000.0
seed  = 1

n   = N + Gamma
Gt  = Gamma / N                    # 0.1
fac = 4 * (Gt + 2) * N**2          # E = fac * lambda = 84000 * lambda
E0  = Gamma**2                     # gap: E_0 = Gamma^2 = 100
print(f"energy rescale E = {fac:.0f} * lambda   gap E0 = {E0}")

# ---- sample the Wishart spectrum and rescale to the edge energy variable ----
rng = np.random.default_rng(seed)
E_all = []
for _ in range(reps):
    M = rng.standard_normal((N, n)) + 1j * rng.standard_normal((N, n))
    w = np.linalg.eigvalsh(M.conj().T @ M).real / (4 * N)     # Fig-4 lambda
    E = fac * w[w > 1e-6 * w.max()]                           # nonzero -> energy
    E_all.append(E[E < 1.05 * Emax])
E_all = np.concatenate(E_all)

# ---- analytic edge densities (mu = hbar = 1) ----
def rho0(E):
    E = np.asarray(E, float); out = np.zeros_like(E); m = E > E0
    out[m] = np.sqrt(E[m] - E0) / (2 * np.pi * E[m]); return out

def rho(E):
    E = np.asarray(E, float); xi = np.sqrt(np.clip(E, 1e-12, None))
    return 0.25 * (jv(Gamma, xi)**2 + jv(Gamma + 1, xi)**2
                   - (2 * Gamma / xi) * jv(Gamma, xi) * jv(Gamma + 1, xi))

# ---- empirical density of states: counts / (reps * binwidth) ----
bins = np.linspace(0, Emax, 220); binw = bins[1] - bins[0]
cent = 0.5 * (bins[:-1] + bins[1:])
dos = np.histogram(E_all, bins=bins)[0] / (reps * binw)

# ---- figure ----
fig, ax = plt.subplots(figsize=(8.6, 5.4))
ax.bar(cent, dos, width=binw, color="#5b8bd0", alpha=0.85, edgecolor="none",
       label=f"Wishart data ({reps} samples)")
Es = np.linspace(0, Emax, 4000)
ax.plot(Es, rho(Es),  color="k",       lw=2.0,               label=r"exact $\rho(E)$ (Bessel, eq. 26)")
ax.plot(Es, rho0(Es), color="#c0392b", lw=1.8, ls=(0, (4, 3)), label=r"leading $\rho_0(E)$ (eq. 37)")
ax.axvline(0, color="#c0392b", lw=3, alpha=0.7)
ax.plot([], [], color="#c0392b", lw=3, label=rf"BPS states ($\Gamma={Gamma}$ at $E=0$)")
ax.axvline(E0, color="gray", ls=":", lw=1.2)
ax.text(E0 + 40, 0.9 * dos.max(), rf"gap $E_0=\Gamma^2={E0}$", color="gray", fontsize=9)
ax.set_xlabel(r"$E$"); ax.set_ylabel(r"$\rho(E)$")
ax.set_xlim(0, Emax); ax.set_ylim(0, None)
ax.set_title(r"Johnson Fig. 5 reproduction: Fig. 4 data zoomed into the hard edge"
             + "\n" + rf"$E=4(\tilde\Gamma+2)N^2\lambda$, the Bessel kernel resolves individual levels",
             fontsize=10.5)
ax.legend(fontsize=9, loc="upper right")
fig.tight_layout()
fig.savefig("johnson_fig5.pdf", bbox_inches="tight")
fig.savefig("johnson_fig5.png", dpi=150, bbox_inches="tight")
print("saved johnson_fig5.pdf / .png")
