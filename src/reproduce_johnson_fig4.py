"""
Reproduction of Figure 4 of C. V. Johnson,
"Fortuitous Chaos, BPS Black Holes, and Random Matrices", arXiv:2601.17122, Sec. VI.

Wishart-type model: M is an N x (N+Gamma) complex Gaussian matrix; H = M^dag M is an
(N+Gamma) x (N+Gamma) positive Hermitian matrix with Gamma exact zero eigenvalues.
The rescaled eigenvalues lambda = lambda_tilde / (4N) follow the Marchenko-Pastur law

    rho_0^MP(lambda) = sqrt((lam_+ - lambda)(lambda - lam_-)) / (pi lambda),
    lam_pm = Gamma_tilde/2 + 1 +/- sqrt(Gamma_tilde + 1),   Gamma_tilde = Gamma/N.

Entry convention: each real and imaginary part is N(0,1) (Johnson's e^{-y^2/2}), i.e.
E|M_ij|^2 = 2. This is fixed by requiring lam_- = 0.001191 at Gamma=10, N=100; then
lam_+ = 2.099 follows automatically.

Usage:  python3 reproduce_johnson_fig4.py
Requires: numpy, matplotlib.  (Johnson uses 1e5 samples; 8e3 already reproduces the curve.)
"""
import numpy as np
import matplotlib.pyplot as plt

# ---- parameters (Johnson Fig. 4) ----
N     = 100     # M has N rows
Gamma = 10      # H = M^dag M has Gamma exact zeros
reps  = 8000    # realizations (Johnson: 1e5; 8e3 is plenty for the leading curve)
seed  = 0

n = N + Gamma                                   # 110: M is N x n, H is n x n
rng = np.random.default_rng(seed)

# ---- sample the Wishart eigenvalues ----
nonzero = []
n_zeros = 0
tol = 1e-6
for _ in range(reps):
    M = rng.standard_normal((N, n)) + 1j * rng.standard_normal((N, n))  # E|M_ij|^2 = 2
    w = np.linalg.eigvalsh(M.conj().T @ M).real                          # n eigenvalues
    w = w / (4 * N)                                                      # lambda = lambda_tilde/(4N)
    thr = tol * w.max()
    nonzero.append(w[w > thr])
    n_zeros += int(np.sum(w <= thr))
lam = np.concatenate(nonzero)
print(f"reps={reps}  nonzero eigenvalues={lam.size}  mean exact zeros/matrix={n_zeros/reps:.2f} (expect {Gamma})")

# ---- Marchenko-Pastur overlay (Johnson eq. 36) ----
Gt = Gamma / N                                  # 0.1
lam_p = Gt / 2 + 1 + np.sqrt(Gt + 1)
lam_m = Gt / 2 + 1 - np.sqrt(Gt + 1)
print(f"Gamma_tilde={Gt}  lam_-={lam_m:.6f}  lam_+={lam_p:.6f}  (Johnson: lam_-=0.001191)")

def rho_mp(x):
    x = np.asarray(x, float)
    out = np.zeros_like(x)
    m = (x > lam_m) & (x < lam_p)
    out[m] = np.sqrt((lam_p - x[m]) * (x[m] - lam_m)) / (np.pi * x[m])
    return out

xx = np.linspace(lam_m, lam_p, 20000)
print(f"integral of rho_0^MP = {np.trapezoid(rho_mp(xx), xx):.4f} (should be 1)")

# ---- figure ----
fig, ax = plt.subplots(figsize=(8.2, 5.4))
ax.hist(lam, bins=140, density=True, range=(0, lam_p * 1.02),
        color="#5b8bd0", alpha=0.85, edgecolor="none",
        label=f"Wishart data ({reps} samples)")
xs = np.linspace(lam_m, lam_p, 800)
ax.plot(xs, rho_mp(xs), color="k", lw=2.2, label="Marchenko--Pastur (eq. 36)")
ax.set_xlabel(r"$\lambda$"); ax.set_ylabel(r"$\rho_0(\lambda)$")
ax.set_xlim(0, lam_p * 1.02); ax.set_ylim(0, None)
ax.set_title(rf"Johnson Fig. 4 reproduction:  $\Gamma={Gamma}$, $N={N}$, "
             rf"$\lambda=\tilde\lambda/(4N)$", fontsize=11)
ax.legend(fontsize=9.5, loc="upper right")

# inset: zoom the origin -> the hard-edge gap. Use TRUE (full-spectrum) density,
# NOT density=True over the sub-range (that would renormalize to the window).
axi = ax.inset_axes([0.12, 0.42, 0.40, 0.48])
ie = np.linspace(0, 0.03, 60); binw = ie[1] - ie[0]
cnt, _ = np.histogram(lam, bins=ie)
axi.bar(ie[:-1], cnt / (lam.size * binw), width=binw, align="edge",
        color="#5b8bd0", alpha=0.85, edgecolor="none")
xi = np.linspace(lam_m, 0.03, 400); axi.plot(xi, rho_mp(xi), color="k", lw=2)
axi.axvline(lam_m, color="#c0392b", ls=":", lw=1.4)
axi.text(lam_m * 1.2, 4.5, rf"gap: $\lambda_-={lam_m:.4f}$", color="#c0392b", fontsize=8)
axi.set_xlim(0, 0.03); axi.set_title("near origin: the hard-edge gap", fontsize=8)
axi.tick_params(labelsize=7)

fig.tight_layout()
fig.savefig("johnson_fig4.pdf", bbox_inches="tight")
fig.savefig("johnson_fig4.png", dpi=150, bbox_inches="tight")
print("saved johnson_fig4.pdf / .png")
