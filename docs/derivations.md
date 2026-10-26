# Derivations

*Living document. Each entry states assumptions, definitions, the derivation, checks, and what is proven vs
conjectured. Status tags: **[derived]** (follows from stated inputs), **[approx]** (uses an uncontrolled
approximation, named), **[num]** (numerical comparison), **[conj]**.*

---

## D1. Decoder variance from the $\mathcal{N}=2$ super-Schwarzian (LMRS zero-energy two-point function)

*2026-09-28. Code: `src/lmrs_predictions.py`; comparison: `scripts/compare_lmrs_variance.py` →
`results/data/lmrs_variance_2026-09-28/`, `results/figures/lmrs_variance_2026-09-28/fig_r_vs_lmrs.pdf`.
Answers the analytic half of Q1 in `docs/research_questions.md`.*

### What is imported and what is derived here

| Ingredient | Source |
|---|---|
| Zero-energy wavefunctions $\lvert Z_j\rangle$, their measure; neutral operator acts as $e^{-\Delta\ell}$ | **Imported:** LMRS eqs. 22–23, 54 (super-Liouville quantization, not redone here) |
| $\mathrm{Tr}\,P_j = e^{S_0}\cos\pi j$; TFD structure at $\beta\to\infty$ | **Imported:** LMRS eqs. 37–39, 68 |
| SYK↔Schwarzian normalization $b_{\hat q}^2(2\alpha_S N)^{-2\Delta}$, value of $\alpha_S$ | **Imported:** LMRS eqs. 63–64, 78–80 |
| Our variance = zero-energy two-point function (Step 1) | **Derived** (exact) |
| Single R-charge sector, $j(N')$, normalization by $\mathrm{Tr}\,P_j$ (Step 2) | **Derived**, verified by exact BPS counts |
| Same-site $n_i$ ≈ hopping bilinear at leading order (Step 3) | **Argued** (Wick); tested numerically, 6–14 % difference |
| Reduction to $\langle e^{-\Delta\ell}\rangle_j$; $\Delta$- and $j$-dependence in closed form (Step 4a–d) | **Derived** analytically from eq. 22; reproduces LMRS eq. 66 |
| Why the prediction overshoots at $N'\le16$ | **Derived** (time-scale analysis) plus toy model; see below |
| LMRS eq. 77 normalization ($2J/N$ should be $2J/N^2$) | Found here, by comparison with FGMS eq. 5.4 |

### Definitions

- One-flavor $\mathcal{N}=2$ SYK with $N'$ complex fermions and supercharge
  $Q = \sum_{i<j<k} C_{ijk}\,\psi^i\psi^j\psi^k$ ($\hat q = 3$), $H = \{Q,\bar Q\}$. This is LMRS eq. (76)
  (the Fu–Gaiotto–Maldacena–Sachdev model).
  - Our code builds $Q$ as the wedge $\omega\wedge$ on $\Lambda^P(\mathbb{C}^{N'})$ (`src/q_scan.form_wedge`).
  - That is the same operator up to charge conjugation $\psi\leftrightarrow\bar\psi$, which maps $j\to-j$. All
    formulas below are even in $j$.
- $P = \lfloor N'/2\rfloor$ (fermion number). $B$ is the BPS space of the $N'$-fermion theory at fermion number $P$,
  with $d = \dim B$. $P_B$ is its projector, and $\tau = \mathrm{Tr}_B/d$.
- Decoder $X = P_B\,\Pi\,P_B$, with $\Pi = 1 - n_i$ and $n_i = \bar\psi_i\psi_i$. Its eigenvalues are
  $\lambda\in[0,1]$.
- UV mean $b = \mathrm{Tr}_D\Pi/D = 1 - P/N'$, and
  $$r \equiv \frac{\mathrm{Var}(\lambda)}{b(1-b)} \le 1 \qquad (\text{since } \lambda\in[0,1]).$$
- R-charge: $\hat q\, j = P - N'/2$ ($\psi$ carries $1/\hat q$).
- LMRS inputs:
  - $\Delta$ is the conformal dimension.
  - $b_{\hat q} = \left[\tan(\pi/2\hat q)/(2\pi)\right]^{1/\hat q}$ (eq. 78).
  - The Schwarzian coupling is $C = \alpha_S N/J$, with $\alpha_S = 0.00842$ at $\hat q = 3$ (eqs. 79–80).

### Step 1: our variance is a zero-energy two-point function [derived, exact]

Let $\mu = \tau(P_B\,n\,P_B)$ and $O = n_i - \mu$. Since $\Pi = 1-n_i$, $\mathrm{Var}(\lambda)$ is the same for
$P_B\Pi P_B$ and $P_B n_i P_B$, and

$$
\mathrm{Var}(\lambda) = \tau\big((PnP)^2\big) - \tau(PnP)^2 = \tau(P\,O\,P\,O\,P)
= \frac{\mathrm{Tr}\,[P_B\,O\,P_B\,O]}{\mathrm{Tr}\,P_B}.
$$

Because $P_B = \lim_{\beta\to\infty} e^{-\beta H}$ restricted to charge $P$, this is LMRS's long-time correlator
$\langle\hat O\hat O\rangle$ with $\hat O = POP$ (their eq. 1), normalized by the number of BPS states in the sector.

### Step 2: $B$ is a single R-charge sector [derived + verified]

Fermion number is conserved and fixes $j$. So $B = \operatorname{ran} P_j$ with $j = (P - N'/2)/\hat q$:

- **$N'$ even:** $P = N'/2$, so $j = 0$.
- **$N'$ odd:** $P = (N'-1)/2$, so $j = -1/(2\hat q) = -1/6$ (LMRS footnote 6: odd $N$ shifts R-charges by
  $1/(2\hat q)$).

**Check** (exact counts, `q_scan.harmonic`, seed 0):

| $N'$ | BPS dimensions by $j$ |
|---|---|
| 8 | $j=0$: 54; $j=\pm\tfrac13$: 27, 27 |
| 10 | $j=0$: 162; $j=\pm\tfrac13$: 81, 81 |
| 12 | $j=0$: 486; $j=\pm\tfrac13$: 243, 243 |
| 9 | $j=\pm\tfrac16$: 81, 81 (plus 3 at $\lvert j\rvert = \tfrac12$) |
| 11 | $j=\pm\tfrac16$: 243, 243 |

- For even $N'$, the ratio $d_{\pm 1/3}/d_0 = 1/2 = \cos(\pi/3)/\cos 0$ is exact, as LMRS eq. (82)
  ($\mathrm{Tr}\,P_j \propto \cos\pi j$) predicts.
- The 3 states at $\lvert j\rvert = 1/2$ for $N'=9$ lie on the boundary the Schwarzian excludes. This is a
  finite-size effect; they are absent at $N'=11$.
- LMRS sum over all $j$ with weight $\cos\pi j/\hat L$. We need only one sector, so we divide by
  $\mathrm{Tr}\,P_j$, not by $N_{\rm BPS}$.

### Step 3: the operator [approx: large $N$, conformal matching]

LMRS eq. (85) treats the neutral bilinear $\psi^k\bar\psi^i$ ($i\neq k$) as a single matter field of dimension
$\Delta = 2\cdot\frac{1}{2\hat q} = \frac13$, with conformal normalization $b_{\hat q}^2$ (two fermion propagators).
For $i = k$,

$$
\langle n_i(t)\,n_i(0)\rangle - \langle n_i\rangle^2 = \langle\bar\psi_i(t)\psi_i(0)\rangle\,\langle\psi_i(t)\bar\psi_i(0)\rangle + O(1/N).
$$

This is the single cross-contraction, identical at leading order to the $i\neq k$ case. So $O$ has the same
conformal two-point function, $b_{\hat q}^2/(Jt)^{2/3}$.

At fixed charge, $\sum_i O_i = 0$ exactly. This affects only the $O(1/N)$ off-diagonal covariance
$c = -1/(N'-1)$, not the leading diagonal term.

### Step 4: the zero-energy average, derived from the LMRS wavefunctions

Rather than quoting LMRS eq. (66), we derive its $\Delta$- and $j$-dependence from their zero-energy wavefunctions.
This shows which inputs the result actually rests on.

**4a. Reduction to a normalized expectation value [derived from the structure of LMRS §2.3–2.4].**
- In sector $j$, the $\beta\to\infty$ limit of the two-sided (TFD) state keeps only its zero-energy component
  $A_{z,j}\,|Z_j\rangle$ (LMRS eqs. 37, 39). Continuum states carry $e^{-E_{s,j}u}$ and drop out as $u,u'\to\infty$.
- A neutral matter operator acts on the Liouville wavefunctions by multiplication by $e^{-\Delta\ell}$ [imported:
  LMRS's bulk treatment, eq. 54]. Hence
  $\mathrm{Tr}[P_j\hat O P_j\hat O] = |A_{z,j}|^2\,\langle Z_j|e^{-\Delta\ell}|Z_j\rangle$ (the structure of eq. 59).
- At $\Delta = 0$ this must equal $\mathrm{Tr}\,P_j$ (eq. 68), so $|A_{z,j}|^2\langle Z_j|Z_j\rangle = \mathrm{Tr}\,P_j$.
  Dividing:
$$
\frac{\mathrm{Tr}[P_j\hat O P_j\hat O]}{\mathrm{Tr}\,P_j} = \frac{\langle Z_j|e^{-\Delta\ell}|Z_j\rangle}{\langle Z_j|Z_j\rangle}.
$$
- Neither $e^{S_0}$ nor the TFD coefficient enters. Only the *shape* of the zero-energy wavefunction matters.

**4b. The wavefunction density [imported: LMRS eq. 22, measure eq. 23].**
- $|Z_j\rangle$ has two fermionic components, $g_1 \propto e^{-\ell/2}K_{\frac12-j}(2e^{-\ell/2})$ and
  $g_2 \propto e^{-\ell/2}K_{\frac12+j}(2e^{-\ell/2})$, with unit-modulus phases and equal prefactors.
- For an operator that depends only on $\ell$, the $a$-integral in the measure (23) is trivial, so
$$
\langle Z_j|f(\ell)|Z_j\rangle \propto \int_{-\infty}^{\infty}\! d\ell\; f(\ell)\; e^{-\ell}\Big[K_{\frac12-j}^2 + K_{\frac12+j}^2\Big]\!\big(2e^{-\ell/2}\big).
$$

**4c. The Mellin integral [derived].** Substitute $x = 2e^{-\ell/2}$, so $d\ell = -2\,dx/x$, $e^{-\ell} = x^2/4$ and
$e^{-\Delta\ell} = (x/2)^{2\Delta}$. Then
$$
I_\nu(\Delta) \equiv \int d\ell\; e^{-(1+\Delta)\ell}\, K_\nu\big(2e^{-\ell/2}\big)^2 = 2^{-1-2\Delta}\int_0^\infty x^{1+2\Delta}K_\nu(x)^2\,dx .
$$
Use the Mellin transform of $K_\nu^2$ (Gradshteyn–Ryzhik 6.576.4 with equal arguments and orders; checked
numerically to $\le 10^{-11}$ for six $(s,\nu)$ pairs):
$$
\int_0^\infty x^{s-1}K_\nu(x)^2\,dx = \frac{\sqrt\pi\,\Gamma(\tfrac s2)\,\Gamma(\tfrac s2-\nu)\,\Gamma(\tfrac s2+\nu)}{4\,\Gamma(\tfrac{s+1}{2})},
\qquad \operatorname{Re}s > 2|\nu| .
$$
With $s = 2+2\Delta$:
$$
I_\nu(\Delta) = 2^{-1-2\Delta}\,\frac{\sqrt\pi}{4}\,\frac{\Gamma(1+\Delta)\,\Gamma(1+\Delta-\nu)\,\Gamma(1+\Delta+\nu)}{\Gamma(\tfrac32+\Delta)} .
$$
Sum the two components, $\nu = \tfrac12\mp j$:
$$
\Gamma(\Delta+\tfrac12+j)\Gamma(\Delta+\tfrac32-j) + \Gamma(\Delta+\tfrac12-j)\Gamma(\Delta+\tfrac32+j)
= (2\Delta+1)\,\Gamma(\Delta+\tfrac12+j)\,\Gamma(\Delta+\tfrac12-j),
$$
using $\Gamma(z+1) = z\Gamma(z)$ on each second factor. Next, $(2\Delta+1)/\Gamma(\Delta+\tfrac32) = 2/\Gamma(\Delta+\tfrac12)$,
and the Legendre duplication formula $\Gamma(\Delta)\Gamma(\Delta+\tfrac12) = 2^{1-2\Delta}\sqrt\pi\,\Gamma(2\Delta)$
gives $2^{-2\Delta}\sqrt\pi/\Gamma(\Delta+\tfrac12) = \Gamma(\Delta)/(2\Gamma(2\Delta))$. Together with
$\Gamma(1+\Delta) = \Delta\Gamma(\Delta)$:
$$
\mathcal N_j(\Delta) \equiv I_{\frac12-j}(\Delta) + I_{\frac12+j}(\Delta)
= \frac18\,\frac{\Delta\,\Gamma(\Delta)^2\,\Gamma(\Delta+\tfrac12+j)\,\Gamma(\Delta+\tfrac12-j)}{\Gamma(2\Delta)} .
$$
The convergence condition $\operatorname{Re}s > 2\nu$ at $\Delta = 0$ is $|j| < \tfrac12$, which is exactly LMRS's
normalizability condition for $|Z_j\rangle$.

**4d. Normalize [derived].** As $\Delta\to0$, $\Delta\Gamma(\Delta)\to1$ and $\Gamma(\Delta)/\Gamma(2\Delta)\to2$, so
$$
\mathcal N_j(0) = \tfrac14\,\Gamma(\tfrac12+j)\,\Gamma(\tfrac12-j) = \frac{\pi}{4\cos\pi j} .
$$
This reproduces the $1/\cos\pi j$ of LMRS's norm, eq. (26), an independent check of step 4b. Therefore
$$
\langle e^{-\Delta\ell}\rangle_j = \frac{\mathcal N_j(\Delta)}{\mathcal N_j(0)}
= \frac{\cos\pi j}{2\pi}\;\frac{\Delta\,\Gamma(\Delta)^2\,\Gamma(\Delta+\tfrac12+j)\,\Gamma(\Delta+\tfrac12-j)}{\Gamma(2\Delta)} ,
$$
which is exactly LMRS eq. (66) divided by $\mathrm{Tr}\,P_j = e^{S_0}\cos\pi j$ (eqs. 38, 68). It is derived here,
not assumed.

**4e. SYK normalization [imported: LMRS eqs. 63–64, 78, 79].**
- In Schwarzian units the operator is normalized to $1/u^{2\Delta}$ at short times, with $u$ measured in units of
  $2C$.
- The SYK conformal correlator is $b_{\hat q}^2/(Jt)^{2\Delta}$.
- So the SYK operator equals $b_{\hat q}^2(2CJ)^{-2\Delta}$ times the Schwarzian one, and $2CJ = 2\alpha_S N'$.
  This is the same step LMRS use in eq. (85).

Combining 4a–4e:

$$
\boxed{\;\mathrm{Var}_j(N') = b_{\hat q}^2\,(2\alpha_S N')^{-2\Delta}\;\frac{\cos\pi j}{2\pi}\;
\frac{\Delta\,\Gamma(\Delta)^2\,\Gamma(\Delta+\tfrac12+j)\,\Gamma(\Delta+\tfrac12-j)}{\Gamma(2\Delta)},
\qquad \Delta = \tfrac13\;}
$$

$$
r_{\rm LMRS}(N') = \frac{\mathrm{Var}_j}{b(1-b)},\qquad
b(1-b) = \begin{cases} \tfrac14 & N' \text{ even}\\[2pt] \dfrac{N'^2-1}{4N'^2} & N' \text{ odd.}\end{cases}
$$

Numerically, at $\hat q = 3$:

$$
\mathrm{Var}_0 = 0.07298\,(0.01684\,N')^{-2/3},\qquad
\frac{\mathrm{Var}_{\pm 1/6}}{\mathrm{Var}_0} = \frac{\cos(\pi/6)\,\Gamma(1)\,\Gamma(2/3)}{\Gamma(5/6)^2} = 0.9204 .
$$

### Checks

1. **Transcription.** The same code reproduces the LMRS Table 1 Schwarzian column:
   - $\psi$ at $j=0$: 0.1112 (table 0.111).
   - $\bar\psi\psi$ at $j=0$: 0.0874 (table 0.0874).
2. **$\Delta\to0$.** The $j$-factor tends to 1 (numerically $1 - 4\times10^{-8}$ at $\Delta = 10^{-8}$, for $j = 0$
   and $\pm\frac16$). The identity operator has $\tau_j(\mathbb 1) = 1$, so dividing eq. (66) by
   $\mathrm{Tr}\,P_j$ is correctly normalized.
3. **Coupling independence.** $J$ cancels ($2CJ = 2\alpha_S N$). This matches the numerics, where $P_B$ does not
   depend on the overall scale of the couplings.
4. **Bound $r\le1$.** The prediction violates it for $N'\le 8$ ($r_{\rm LMRS} = 1.11$ at $N'=8$). There the
   Schwarzian regime (gap $\propto 1/C \propto 1/N$) is not reached, so the formula should only be trusted at
   larger $N'$.
5. **Parity.** At fixed $N'$, the prediction makes odd $N'$ less rigid (factor 0.92). This is the observed even/odd
   alternation, e.g. $r(13) = 0.617 < r(14) = 0.642$.

### Comparison with exact numerics [num]

Data: $\hat q = 3$, campaign `summary.csv`, $N' = 5\ldots15$, with 60 down to 3 realizations per point.

| $N'$ | 9 | 10 | 11 | 12 | 13 | 14 | 15 |
|---|---|---|---|---|---|---|---|
| $r$ (SYK) | 0.728 | 0.746 | 0.666 | 0.687 | 0.617 | 0.642 | 0.577 |
| $r_{\rm LMRS}$ (no free parameter) | 0.957 | 0.957 | 0.833 | 0.847 | 0.744 | 0.765 | 0.675 |
| ratio | 0.761 | 0.780 | 0.799 | 0.811 | 0.830 | 0.840 | 0.854 |
| Haar: $r = a$ | 0.643 | 0.643 | 0.526 | 0.526 | 0.425 | 0.425 | 0.340 |

- **The parameter-free prediction does not match as it stands.** It overshoots by 15–24 % (rms 0.17). The
  discrepancy shrinks monotonically with $N'$.
- **The ratio is smooth in $N'$.** The even/odd alternation of $r$ is removed entirely by the $j$-dependence, with
  no parameters. Without it, the odd points would sit 8 % lower. This is the one sharp parameter-free success.
- **One-parameter fit (a fitted ansatz, not derived).** $r/r_{\rm LMRS} = 1 - c/N'$ with $c = 2.20$ fits
  $N' = 9\ldots15$ with max residual 0.006. A two-parameter fit $A - B/N'$ gives $A = 0.974$, $B = 1.93$.
  - A $1/N'$ correction is generically expected: the relevant times are $t\sim C\sim N/J$, and $\alpha_S$ itself has
    $1/N$ corrections.
  - But $c = 2.2$ (a 15 % effect at $N' = 15$) is unexplained, so the corrected curve is not a prediction.
- **Does the data prefer the Schwarzian $N'^{-2/3}$ shape? No — not at these sizes.** Fits over $N' = 9\ldots15$,
  weighted by the s.e. (`scripts/compare_lmrs_variance.py` → `fit.json`, `alternative_models`):

  | model for $r(N')$ | free params | $\chi^2/\text{dof}$ |
  |---|---|---|
  | $r_{\rm LMRS}$ | 0 | — (rms 0.17) |
  | $A\, r_{\rm LMRS}$ (equivalently, $\alpha_S$ free) | 1 | 335 |
  | $r_{\rm LMRS}(1 - c/N')$, $c = 2.20$ | 1 | 9.0 |
  | $K N'^{-p}\times$ (LMRS parity factor), $p = 0.44$ | 2 | 1.3 |
  | $K N'^{-p}(1 - e\,[N'\text{ odd}])$, $p = 0.45$, $e = 0.069$ | 3 | 0.2 |
  | $K N'^{-p}$, no parity | 2 | 700 |
  | Haar $a(1 + c/N')$ | 1 | 3963 |

  - The pure $N'^{-2/3}$ shape fails even with its normalization free ($\chi^2/\text{dof} = 335$).
  - "$N'^{-2/3}$ with a $1/N'$ correction" and "$N'^{-0.44}$" describe the data equally well. Seven points cannot
    separate an asymptotic power from a large correction.
  - The fitted parity amplitude (0.069) is close to the LMRS $j$-factor value (≈ 0.076), which supports the R-charge
    structure independently of the $N'$-dependence.
  - Every Haar-based model fails badly.
- **$\alpha_S$ sensitivity.** $r \propto \alpha_S^{-2/3}$, so $A = 0.974$ corresponds to $\alpha_S \approx 0.0088$.
  That is between LMRS's numerical 0.00842 and the large-$\hat q$ extrapolation 0.0092 (LMRS footnote 5). With
  $\alpha_S = 0.0092$ the one-parameter fit degrades (max residual 0.02).
- **Local slopes.** $d\ln r/d\ln N' = -0.41\ldots-0.47$. That is not $-2/3$, but it is what
  $r_{\rm LMRS}(1 - 2.2/N')$ implies at $N'\approx12$:
  $$-\tfrac23 + \frac{c/N'}{1 - c/N'} = -0.44 .$$
- **Out-of-sample prediction** (made before the run finished): $r(N'=16) = 0.597$–$0.604$ (one-/two-parameter
  $1/N'$ fits).
  - Caveat: naive extrapolation of the even-$N'$ data also gives $\approx 0.61$, so this test is weak.
- **What LMRS did.** They compare at a single size, $N = 16$, using the large-$N$ $\alpha_S$, and apply no finite-$N$
  correction. They describe the result as "agreement, within a few percent" (introduction). Their neutral
  $\bar\psi_i\psi_j$ ($i\ne k$, all $j$) is 0.079 ± 0.001 measured vs 0.0874 predicted: 10 % low, outside the
  statistical error, and not discussed. At $N' = 16$ our same-mode operator in the $j=0$ sector is expected to be
  about 14 % low. So at the one size LMRS tested, our discrepancy is the same order as theirs. LMRS give no
  $N$-dependence, so they cannot say whether the approach is $1/N$-like.

### Out-of-sample test: R-charge dependence at fixed $N'$ [num, pre-registered]

*`scripts/lmrs_jsector_test.py` → `results/data/lmrs_variance_2026-09-28/jsector.jsonl`; seeds 0–3, 4 modes,
peak 0.8 GB.*

For even $N'$, BPS states also exist at $P = N'/2 \pm 1$ ($j = \pm\tfrac13$). With the same couplings, Step 4
predicts, with no free parameter and before any measurement:
$$
\frac{\mathrm{Var}_{1/3}}{\mathrm{Var}_0} = \frac{\cos(\pi/3)\,\Gamma(\tfrac76)\,\Gamma(\tfrac12)}{\Gamma(\tfrac56)^2} = 0.645 \quad(\text{LMRS}),
\qquad \approx 0.66 \quad(\text{UV-capped toy}).
$$
Comparing at fixed $N'$ removes most of the $N'$-dependent finite-size effects that spoil the absolute comparison.

| $N'$ | measured (mean ± s.e.) | LMRS | capped toy | Haar null | no $j$-dependence |
|---|---|---|---|---|---|
| 8 | 0.642 ± 0.009 | 0.645 | 0.663 | 0.586 | 1 |
| 10 | 0.668 ± 0.006 | 0.645 | 0.664 | 0.576 | 1 |
| 12 | 0.677 ± 0.006 | 0.645 | 0.662 | 0.567 | 1 |
| 14 | 0.683 ± 0.004 | 0.645 | 0.660 | 0.560 | 1 |

(The $j = +\tfrac13$ and $-\tfrac13$ sectors agree to all digits. That follows from particle–hole symmetry of the
model, so it is a check of the code, not of LMRS.)

- **Supports LMRS.** A 32 % sector dependence is measured; LMRS predict 35.5 % from Γ-functions alone. The Haar null
  (44 %, moving the other way with $N'$) and "no $j$-dependence" (0 %) are both clearly excluded. The prediction is
  within 6 % at every $N'$.
- **Does not converge.** The measured ratio drifts *away* from 0.645 as $N'$ grows (0.642 → 0.683; 9σ from 0.645 at
  $N'=14$). So finite-$N'$ effects of order 5 % persist in the ratio and do not visibly shrink by $N'=14$. That is
  consistent with the no-conformal-window picture below, where the relevant expansion parameter is $1/(2\alpha_S N')$,
  not $1/N'$. But this test cannot show convergence to the Schwarzian value.

### Why the numbers don't match: no conformal window at $N'\lesssim 16$

*2026-09-28. Code: `src/lmrs_predictions.py` (`zero_energy_expectation`, `var_capped`),
`scripts/lmrs_mismatch_checks.py` → `results/data/lmrs_variance_2026-09-28/mismatch_checks.jsonl`
($N' = 8\ldots14$, seeds 0–3, peak 0.8 GB under memwatch).*

**1. What the zero-energy correlator actually samples [derived].** The LMRS result is an expectation over the
zero-energy wavefunction (eq. 22) in the renormalized length $\ell$:
$$
\frac{\mathrm{Tr}[P_j O P_j O]}{\mathrm{Tr}\,P_j} = \big\langle\, G_{\rm matter}(t(\ell))\,\big\rangle_{Z_j},\qquad
\rho_j(\ell) \propto e^{-\ell}\Big[K_{\frac12-j}^2 + K_{\frac12+j}^2\Big]\!\big(2e^{-\ell/2}\big),\qquad
J t(\ell) = 2\alpha_S N'\, e^{\ell/2}.
$$
With the conformal input $G_{\rm matter} = b_{\hat q}^2 (Jt)^{-2\Delta}$, this integral reproduces the closed form of
Step 4 to $10^{-6}$ for all tested $\Delta$ and $j$ (test `test_wavefunction_integral_reproduces_gamma_formula`).
The physical times it samples are set by $2\alpha_S N'$, and $\alpha_S = 0.00842$ is small:

| $N'$ | median $Jt$ | weight where conformal $> $ UV maximum $\tfrac14$ |
|---|---|---|
| 10 | 0.50 | 67 % |
| 15 | 0.75 | 49 % |
| 20 | 1.00 | 35 % |
| 50 | 2.5 | 4 % |
| 120 | 6.0 | 0 % |

(Weight $\rho_0(\ell)\,e^{-\ell/3}$, $j=0$.) At our sizes, the "zero-energy" correlator is dominated by times
**shorter than the microscopic time $1/J$**, where the conformal propagator is not valid. A clean conformal window
(median $Jt\approx10$) needs $N'\approx 200$. The same fact shows up elsewhere:
- $r_{\rm LMRS} > 1$ for $N'\le 8$.
- The measured BPS gap is $0.19\,J$ at $N'=14$ (not $\ll J$), and $E_{\rm gap}/N'$ drifts ($3.7\to2.6$ in our units),
  where the Schwarzian predicts a constant.

**2. Consequence: the conformal input overshoots [derived sign; size model-dependent].** The true connected
correlator is bounded by its equal-time value, $b(1-b) = \tfrac14$ for $n_i$. It is expected to lie below the
conformal extrapolation at $Jt\lesssim1$ (the leading non-conformal correction in SYK reduces $G$). Replacing the
conformal input by conformal-until-saturation, $\min(G_{\rm conf}, V_{\rm UV})$, gives a toy with **no fitted
parameter** (`var_capped`):

| $N'$ | 9 | 10 | 11 | 12 | 13 | 14 | 15 |
|---|---|---|---|---|---|---|---|
| LMRS (conformal) | 0.236 | 0.239 | 0.207 | 0.212 | 0.185 | 0.191 | 0.168 |
| capped toy | 0.172 | 0.180 | 0.162 | 0.170 | 0.153 | 0.161 | 0.145 |
| measured $\mathrm{Var}$, $n_i$ | 0.180 | 0.187 | 0.165 | 0.172 | 0.153 | 0.161 | 0.144 |
| measured, $c_k^\dagger c_i$ (4 seeds) | 0.203 | 0.207 | 0.182 | 0.189 | 0.164 | 0.170 | — |
| capped toy, $c_k^\dagger c_i$ (cap $\frac{P(N'-P)}{N'(N'-1)}$) | 0.183 | 0.190 | 0.169 | 0.177 | 0.158 | 0.165 | — |

- The toy captures the size of the overshoot, its decrease with $N'$, and the parity structure.
- Its near-exact agreement for $n_i$ at $N'\ge13$ is **partly coincidental**:
  - For the hopping operator it undershoots by 10 % → 3 %.
  - A smooth interpolation $1/(1/G_{\rm conf}+1/V_{\rm UV})$ undershoots by about 40 %.
  - So the result depends on the crossover shape, which the toy does not know.
- Robust statement: for both operators the data lie **between** the conformal LMRS value and the hard-capped
  value, converging toward the latter as $N'$ grows.

**3. Secondary: same-site vs. hopping bilinear [num].** LMRS eq. (85) treats $c_k^\dagger c_i$ ($i\ne k$); we
need $n_i$. They agree at leading order (Step 3), and the measured ratio
$\mathrm{Var}_{\rm hop}/\mathrm{Var}_{n_i} = 1.14,\,1.13,\,1.11,\,1.10,\,1.09,\,1.08,\,1.06$ ($N' = 8\ldots14$) indeed
tends to 1, like $1/N$. At $N'=14$ this accounts for about 6 of the 16 percentage points. The hopping operator is
11 % below LMRS at $N'=14$, consistent with LMRS's own 10 % miss at $N=16$ for the same operator (Table 1).

**4. Ruled out / not the cause.**
- *$\alpha_S$ alone.* A constant rescaling cannot fit ($\chi^2/\text{dof} = 335$). The gap cannot calibrate
  $\alpha_S$ at these $N'$ (item 1).
- *The sum rule $\sum_i O_i = 0$.* It fixes the $O(1/N)$ covariance exactly; its effect on the diagonal variance is
  $O(1/N^2)$, because the zero-energy charge fluctuations across $j$ are $O(1)$.

**Normalization note.** LMRS eq. (77) prints $\langle C\bar C\rangle = 2J/N$; FGMS eq. (5.4), whose conventions
LMRS say they use, has $2J/N^2$, and only $N^2$ gives an extensive $H$. So our couplings ($E|C|^2 = 2$) correspond
to $J = N'^2$. $J$ cancels in the variance prediction; it matters only for the gap check.

**Bottom line [interp].** The super-Schwarzian formula is the leading term of an expansion whose effective small
parameter is $1/(2\alpha_S N')\approx 60/N'$ in time units. At $N'\le16$ that parameter is not small, so the formula
is evaluated outside its regime, and the 15–24 % overshoot is the expected sign and size. The fitted $1 - 2.2/N'$
and the effective exponent $0.44$ (instead of $2/3$) are finite-$N$ crossover artifacts, not evidence against
LMRS. They are not evidence for $N'^{-2/3}$ either: that regime sets in only at $N'\sim 10^2$. A controlled
prediction at our sizes needs the exact large-$N$ (Schwinger–Dyson) matter correlator inserted into the
zero-energy average, taking care not to double count the $h=2$ mode already in the Schwarzian.
**Update: see D2.** The double-scaled super-chord result of Boruch–Lin–Yan, which is UV-complete at finite $\lambda$,
reproduces $r$ to 2 % with no free parameters.
The same caveat applies to any gravitational prediction of $X$ or $e_3$ compared against $N'\le16$ data.

### Multi-mode slots $\Pi_T$, $|T| = k$ [approx, exploratory]

At leading large $N$, different sites factorize, and
$\langle(1-n_i)(t)\,(1-n_i)(0)\rangle = b^2 + g(t)$, with $g$ the bilinear correlator. Hence the connected
correlator of $\Pi_T$ is

$$
\sum_{m=1}^{k}\binom{k}{m}\, b^{2(k-m)}\, g(t)^m .
$$

The term $g^m$ is a $2m$-fermion composite of dimension $m/3$. Treating it as one field (the same approximation as
Step 3):

$$
\mathrm{Var}_k = \sum_{m=1}^{k}\binom{k}{m}\, b^{2(k-m)}\, b_{\hat q}^{2m}\,(2\alpha_S N')^{-2m/3}\, F\!\left(\tfrac m3, j\right),
\qquad b_k = \frac{\binom{N'-k}{P}}{\binom{N'}{P}},
$$

where $F(\Delta, j)$ is the $j$-factor of Step 4.

The comparison is **not a quantitative test**:
- The measured/predicted ratio rises toward 1 with $N'$.
- But $N'(1-\text{ratio}) \approx 3.9$ ($k=2$) and $\approx 6.5$ ($k=3$), versus 2.2 ($k=1$), and a parity wobble
  remains.
- Heavier composites are further from the Schwarzian regime; the predicted $r_k$ even exceeds 1 at $N'\le12$.
- The obvious corrections are exact finite-$N$ combinatorics (sampling without replacement) and the composite
  approximation.

### Status

- **Derived:** Steps 1, 2, and 4a–4d, the last analytically from the LMRS wavefunctions (eq. 22) via a Mellin
  integral and the duplication formula; checks 1–3. **Imported without re-derivation:** the super-Liouville
  quantization behind eq. 22, the TFD structure, and the SYK↔Schwarzian normalization including $\alpha_S$.
- **Approximation:**
  - Step 3, the bilinear as a single $\Delta = 1/3$ field. This is LMRS's own approximation, about 10 % accurate at
    $N=16$ in their Table 1.
  - The value of $\alpha_S$ (numerical, large $N$).
- **Numerical observations:**
  - The parameter-free prediction overshoots by 15–24 %, decreasing with $N'$.
  - The parameter-free R-charge (parity) factor is confirmed.
  - $r_{\rm SYK}/r_{\rm LMRS} = 1 - 2.2/N'$ fits $N' = 9\ldots15$, but $c$ is **fitted, not derived**. Candidate
    sources are the $1/N$ correction to $\alpha_S$, non-conformal corrections at $t\sim N/J$, $O(1/N)$ Wick terms in
    the $i=k$ correlator, and fixed-charge constraints.
- **Consequence for Q1 [conj — partially supported]:**
  - *Established at these sizes:* $r$ is far from the Haar/Wachter value $r = a \propto (\sqrt3/2)^{N'}$ (70 %
    above it at $N'=15$), and it decays slowly, like a power.
  - *Consistent with, not established:* the Schwarzian law $r\propto N'^{-2/3}$. The data's effective exponent is
    0.44, and it matches LMRS only after a fitted $1/N'$ correction.
  - *Sharpest parameter-free agreement:* the even/odd (R-charge $j$) structure.
  - Deciding the asymptotic power needs $N'\gtrsim 20$ (e.g. stochastic trace estimation, `src/q_scan_mf.py`), or
    an analytic computation of the $1/N$ correction.
- **Not computed:** the crossing weight $X$ (zero-energy OTOC/TOC at $\Delta = 1/3$).
  - LMRS give only two limits: $\Delta\gg1$, where the ratio is $\sim e^{-3.52\Delta}$ (eq. 159); and
    $\Delta\to0$, where all contractions have equal weight, i.e. $X\to1$ (eq. 169).
  - The finite-$\Delta$ value needs the propagator integral (LMRS eq. 145). This is the next analytic target, and
    it would test $e_3$ in the chord-invariant note.

---

## D2. Decoder variance and BPS fraction at finite λ (double-scaled super-chords, Boruch–Lin–Yan)

*2026-09-28. Source: Boruch, Lin & Yan, arXiv:2308.16283 ("BLY"; note `literature/notes/boruch_lin_yan_2023.md`).
Code: `src/dssyk_n2.py` (tests `tests/test_dssyk_n2.py`); comparison `scripts/compare_dssyk.py` →
`results/data/dssyk_2026-09-28/comparison.csv`, `results/figures/dssyk_2026-09-28/fig_r_a_vs_chords.pdf`.*

### What is imported and what is done here

| Ingredient | Source |
|---|---|
| Super-chord Hilbert space, HH state, finite-λ two-point function (4.7)/(4.8), BPS fraction (3.17), index (A.4) | **Imported:** BLY (not re-derived) |
| Identification of BLY's normalized $\mathrm{Tr}(\Pi_{0,j}W\Pi_{0,j}W)$ with our $r$ | **Derived** here (below) |
| Parameter map $\lambda = 2p^2/N'$, $\Delta = 1/p$, $j = (P-N'/2)/p$ | **Derived**, from BLY eqs. 1.6, 1.8 and D1 Step 2 |
| (A.4) reproduces all our exact BPS dimensions | **Verified** (test) |
| Comparison with exact $r$, $a$ for $N' = 5\ldots16$ | **Numerical** |

### Identification with $r$ [derived]

- BLY normalize the neutral operator so that its infinite-temperature two-point function in the charge sector is 1
  ($W_LW_R \propto q^{\Delta(n_O+n_X)}$ with unit constant, so at $n=0$, i.e. $|\Omega\rangle$, it equals 1).
- The HH state is $|\Psi,j\rangle \propto \Pi_{0,j}|\Omega,j\rangle$. So
$$
\frac{\langle\Psi,j|W_LW_R|\Psi,j\rangle}{\langle\Psi,j|\Psi,j\rangle}
= \frac{\mathrm{Tr}(\Pi_{0,j}W\Pi_{0,j}W)/\mathrm{Tr}\,\Pi_{0,j}}{\mathrm{Tr}_j(WW)/\mathrm{Tr}_j\mathbb 1} .
$$
- For $W = n_i-\mu$ in $\Lambda^P$, $\mathrm{Tr}_j(WW)/\mathrm{Tr}_j\mathbb 1 = \frac PN(1-\frac PN) = b(1-b)$ exactly,
  and the numerator is $\mathrm{Var}(\lambda)$ (D1 Step 1). Hence
$$
r = \frac{\mathrm{Var}(\lambda)}{b(1-b)} \;\overset{\rm BLY}{=}\;
\frac{(q^{2+4\Delta};q^2)_\infty\,(q^{1\pm2j};q^2)_\infty\,(q^2;q^2)_\infty}{(q^{2+2\Delta};q^2)_\infty^2\,(q^{1+2\Delta\pm2j};q^2)_\infty},
\qquad q = e^{-2p^2/N'},\ \Delta = \tfrac1p .
$$
  (Products over ± are implied.)
- The $\lambda\to0$ limit is $(2\lambda)^{2\Delta}$ times the LMRS $j$-factor of D1 (BLY 4.9; checked numerically,
  ratio $1 - 3\times10^{-3}$ at $\lambda = 0.3$).
- Unlike D1, no conformal normalization $b_{\hat q}$ and no $\alpha_S$ enter. The operator is normalized in the UV,
  which is why this can capture the crossover.

### From $r$ to the decoder second moment $T_2$ [derived; added 2026-10-16]

*Used by the parameter-free $T_4$ prediction (progress report 2026-10-13 §9.3; D7) and by
`scripts/analyze_uplift_heterogeneity.py::chord_T2`.*

- Decoder moments use $\mu=2\lambda-1\in[-1,1]$ (the eigenvalues of $P_BUP_B$, $U=(-1)^{n_i}$), with
  $T_{2m}=\langle\mu^{2m}\rangle$. Expanding,
$$
T_2=4\langle\lambda^2\rangle-4\langle\lambda\rangle+1=4\,\mathrm{Var}(\lambda)+\big(2\langle\lambda\rangle-1\big)^2 .
$$
- **Mean.** On $\Lambda^P$, $\sum_i n_i=P$ exactly. Averaged over the probed site,
  $\frac1{N'}\sum_i\tau(P_B\Pi_iP_B)=1-P/N'=b$ in every realization. Hence $\langle\lambda\rangle=b$ exactly for the
  site-averaged decoder. For a single fixed site (e.g. the last one) it holds in the ensemble mean only.
- **Variance.** $\mathrm{Var}(\lambda)=b(1-b)\,r$ by definition of $r$.
- With $(2b-1)^2=1-4b(1-b)$:
$$
\boxed{\;T_2=1-4b(1-b)(1-r)\;}\qquad\Longleftrightarrow\qquad r=\frac{T_2-c^2}{1-c^2},\quad c=2b-1 .
$$
- Even $N'$: $b=\tfrac12$, so $T_2=r$. Odd $N'$: $c=1/N'$, a small offset.
- **Status.** Exact given $r$ (site-averaged). The only approximation in the parameter-free $T_2$ is the BLY value
  of $r$ above, evaluated by `src/dssyk_n2.r_chord` (= `two_point_chord(2p²/N′, (P−N′/2)/p, 1/p)`) and tabulated as
  `r_chord` in `results/data/dssyk_2026-09-28/comparison.csv`; it matches the exact $r$ to about 2 %.

### Results [num] (no free parameters)

| $N'$ | $j$ | $\lambda$ | $r$ (SYK) | $r$ chord | $r$ LMRS | SYK/chord | $a$ exact | $a$ chord (3.17) |
|---|---|---|---|---|---|---|---|---|
| 6 | 0 | 3.00 | 0.905 | 0.913 | 1.345 | 0.990 | 0.900 | 0.900 |
| 8 | 0 | 2.25 | 0.819 | 0.832 | 1.110 | 0.984 | 0.771 | 0.789 |
| 10 | 0 | 1.80 | 0.746 | 0.756 | 0.957 | 0.987 | 0.643 | 0.671 |
| 12 | 0 | 1.50 | 0.687 | 0.690 | 0.847 | 0.996 | 0.526 | 0.559 |
| 13 | ∓1/6 | 1.38 | 0.617 | 0.613 | 0.744 | 1.007 | 0.425 | 0.456 |
| 14 | 0 | 1.29 | 0.642 | 0.635 | 0.765 | 1.011 | 0.425 | 0.459 |
| 15 | ∓1/6 | 1.20 | 0.577 | 0.565 | 0.675 | 1.021 | 0.340 | 0.371 |
| 16 | 0 | 1.13 | 0.601 ± 0.005 | 0.589 | 0.700 | 1.022 | 0.340 | 0.373 |
| 14 | ±1/3 | 1.29 | 0.449 | 0.418 | 0.504 | 1.074 | 0.243 | 0.265 |

(Full table, all $N' = 5\ldots16$ and $j = \pm\tfrac13$ at $N' = 8\ldots14$: `comparison.csv`.)

- **$r$, main sectors:** agreement within 2.2 % for every $N' = 5\ldots16$, including the parity alternation. The
  super-Schwarzian misses by 15–60 %.
- **Residual:** a systematic drift of SYK/chord from 0.98 ($N'=7$) to 1.02 ($N'=16$). In the $j = \pm\tfrac13$
  sectors it goes from 0.98 to 1.07 ($N' = 8 \to 14$). The fixed-$N'$ ratio $\mathrm{Var}_{1/3}/\mathrm{Var}_0$ is
  0.645 in the chord formula as well (near-coincidence: the chord $r$-ratio varies with $\lambda$ but compensates
  the UV-variance ratio), so the measured drift of that ratio (D1) is **not** explained by the chord formula either.
- **BPS fraction:** (A.4) is exact. (3.17) at $\lambda = 18/N'$ overshoots by up to 10 %, growing with $N'$. This is
  a pure finite-$p$ effect ($\cos^N(\pi/2p)$ vs $e^{-N\pi^2/8p^2}$), and it gives a scale for how large finite-$p$
  corrections to double-scaled formulas are at our sizes.

### Interpretation and caveats

- **[interp]** This confirms D1's diagnosis quantitatively. At $N'\le16$ the relevant boundary times are
  $\lesssim 1/J$, and the super-chord theory, which is the UV completion of the Schwarzian, captures the crossover.
- **[caveat]** The double-scaled theory at $\lambda = 18/N'$ is *not* the $N'\to\infty$ limit at fixed $p = 3$. Its
  Schwarzian coupling is $\alpha_S = 1/(8p^2) = 0.0139$ (BLY: $C = \alpha_S N = 1/4\lambda$), versus 0.00842 for the
  true $\hat q = 3$ model. So the agreement is expected to degrade as $N'$ grows. The positive drift (+2 % at
  $N'=16$) may be its onset. A finite-$p$ correction to the chord rules is the natural next refinement. For
  example, the exact crossing weight of a single-site bilinear with a random 3-body chord is $1 - 2p/N$ instead of
  $q^\Delta = e^{-2p/N}$; untested.
- **[caveat: λ identification]** We use BLY's definition $\lambda = 2p^2/N'$ (eq. 1.6), fixed before comparing. The
  prediction scales as $r\propto\lambda^{0.39\ldots0.58}$ ($N' = 8\ldots16$), so a 5 % change in the finite-$p$
  identification of $\lambda$ shifts $r$ by about 2.5 %. The 2 % agreement, and certainly the $-2\%\to+2\%$ drift, is
  therefore at the resolution of the double-scaled mapping at $p = 3$. Not every "reasonable" finite-$N$ choice is
  equivalent. Defining $\lambda = -\ln E[(-1)^{|I\cap J|}]$ from exact 3-subset overlaps is undefined for
  $N'\le10$ (the average sign is negative) and gives $r = 0.86$ at $N' = 14$. That choice is not the right finite-$N$
  object for $\mathcal N=2$ chords anyway, whose friend/enemy weights also involve $2^s$ factors.
- **[not done]** The zero-temperature OTOC (our crossing weight $X$) and the 6-letter words ($e_3$). BLY §5 gives
  the one-particle wormhole supercharges and $W_LW_R\sim q^{\Delta(n_O+n_X)_{\rm tot}}$ (5.18), but does not
  evaluate them. For our operators, the matter–matter crossing weight must be the site-local one (1 for disjoint
  sites, since those bilinears commute), not the random-operator $q^{\Delta^2}$. This is the concrete form of action
  item 2.

### Implications for the project goals

**1. Free compression is the zero-length-wormhole term.**

*Provenance.*

| Ingredient | Source |
|---|---|
| $D(j) = P_0$: BPS fraction = probability of zero chord number | **BLY** eq. 1.1, derived 3.12–3.15; closed form 3.17 |
| Zero-T two-point function = $\sum_n P_n q^{2\Delta n}$ (norm sum 3.10 with $q^{2\Delta n}$ inserted; term/norm = $P_n$) | **BLY** eqs. 4.5–4.6 and text after 3.11; closed form 4.7–4.8 |
| That quantity equals our $r$ | **Derived here** (BLY normalization + D1 Step 1) |
| Splitting off $n=0$; $\eta_r = \langle q^{2\Delta n}\rangle_{n\ge1}$ | **Derived here** (trivial algebra) |
| $n=0$ bookkeeping (BLY's summand has a singular $(q^2;q^2)_{-1}$ at $n=0$) | **Ours:** set $P_0 = D(j)$ and check against both closed forms |
| "$P_0$ ↔ Haar/Wachter value" | **Interpretation, ours.** BLY do not discuss freeness; the Haar value $\rho_2 = a$ is free compression |

Write $a_\lambda \equiv D(j)$ (double-scaled) and $r_\lambda$ for BLY 4.8. Then, exactly within the double-scaled theory,
$$
r_\lambda = \underbrace{a_\lambda}_{=P_0} + \sum_{n\ge1} P_n\,q^{2\Delta n},
\qquad
\eta_r^{(\lambda)} = \frac{r_\lambda-a_\lambda}{1-a_\lambda} = \big\langle q^{2\Delta n}\big\rangle_{n\ge1}.
$$
Checked from the HH recursion (BLY 3.2–3.11) for $N' = 8, 12, 14, 16$: with $P_0 = a_\lambda$ assigned, both
$\sum P_n = 1$ (norm 3.11) and $\sum P_n q^{2\Delta n} = r_\lambda$ (4.7) hold to $\le10^{-10}$
(`hh_length_distribution`, test `test_free_compression_is_zero_length_term`).

*Reading [interp].* The zero-length wormhole contributes exactly the value that free compression assigns to $r$ at
that BPS fraction. The deviation from freeness is the contribution of nonzero-length wormholes.

*Caveats.*
- **Equality of values, not of structures.** Only the second cumulant is tested.
- **$a_\lambda$ vs exact $a$.** $a_\lambda$ is 0–10 % larger at $p=3$. The exact Haar statement uses $a$; they agree
  only within the double-scaled approximation. This explains the ~0.04 offset of measured $\eta_r$ from the chord
  curve.
- **$r\ge a$ is automatic in the chord theory, not an exact finite-$N$ theorem.** Individual realizations violate it
  where $a\gtrsim0.9$ and $d$ is tiny ($q=3$, $N'=5,6$; $q=5$, $N'=9,10$; 2/300 at $q=3$, $N'=7$; min $r/a = 0.965$).
  Realization averages satisfy it within errors everywhere.

**[conj → falsified in strong form, 2026-09-29]** "The whole free-compression law is the $\ell=0$ truncation of the
super-chord computation." Beyond $m_2$ the decoder is non-free even at $x=0$; see
`research/notes/decoder_moments_2026-09-29.md`. That note also gives the **exact finite-$N$ version** of the identity:
$E[T_2]=\sum_k w_k\hat\chi_k(U)$ with $w_0=a$ exactly (operator-size decomposition of the BPS projector).

**2. Prediction: $\eta_r$ peaks, then decays [chord, λ = 18/N′].**

| $N'$ | 8 | 12 | 16 | 20 | 24 | 30 | 40 | 60 | 100 |
|---|---|---|---|---|---|---|---|---|---|
| $\eta_r$ chord | 0.20 | 0.30 | 0.34 | 0.361 | 0.363 | 0.35 | 0.32 | 0.25 | 0.18 |
| $\eta_r$ measured | 0.21 | 0.34 | 0.40 | | | | | | |

- The measured deceleration (§4 of the chord-invariant note) is the approach to a maximum near $N'\approx 22$, not
  a plateau.
- The measured values lie about 0.04 above the chord curve. That is consistent with the finite-$p$ excess of
  $a_{\rm chord}$ over the exact $a$ (7–10 %).
- Asymptotically $a\sim e^{-\pi^2 N'/72}\to0$ while $r\sim N'^{-2/3}$. So $r/a\to\infty$ (ever further from Wachter),
  yet $r\to0$ (ever further from the commuting UV). $\eta_r$ is not the natural large-$N'$ measure; $r/a$ and $r$ are.
- The same qualitative turnover follows for fixed $p = 3$ from LMRS, since there $r\propto N'^{-2/3}$ too.

**3. Q1 answered structurally.** In the chord theory the decoder variance is power-law small ($\lambda^{2\Delta}$), and
the Haar value is exponentially small. The accessible-$N'$ data are quantitatively in the chord regime (D2).

**4. What remains for the spectral density.** Higher moments $m_k = \tau((P\Pi P)^k)$ are zero-temperature $k$-point
functions of one fixed operator. Two points need care:
- $n_i$ is a projector, not a random matter operator, so its matter chords are not Gaussian by themselves. Writing
  $n_i = \bar\psi_i\psi_i$ with $\psi_i$ a single-fermion chiral chord ($\Delta_O = 1/p$) restores a Wick structure,
  because the infinite-temperature state is Gaussian.
- Same-site fermion chords then cross with sign $-1$, and disjoint-site bilinears commute (weight 1).

This is the multi-particle wormhole machinery of BLY §5, the same computation that gives $X$ and $e_3$.

---

## D3. Operator size per unit of chord number is $p-1$ (the body number of $H$)

*2026-09-30. Prompted by Lin 2022 (`literature/notes/lin_2022.md`, eqs. 53–59). Code:
`scripts/check_size_per_chord.py` → `results/data/size_per_chord_2026-09-30/` (seeds 0–1; $p=3$: $N'=10,12$;
$p=5$: $N'=10,12,14$; peak memory 1.08 GB under memwatch; commit 166252a + working tree). Existing $p=3$ size data:
`results/data/parity_family_2026-09-29/size_q3_N*.jsonl` (seeds 0–3). Resolves, at the level stated below, the open
question "why does one super-chord pair carry size 2 rather than $p$?" of `research/pdfs/progress_2026-09-29.pdf` §8.*

### Definitions
- One-flavor $\mathcal N=2$ SYK on $N=N'$ complex fermions, $Q=\sum_{|I|=p}C_I\psi_I$ with $p$ odd,
  $H=\{Q,Q^\dagger\}$ on the half-filled sector $\Lambda^F$, $F=\lfloor N/2\rfloor$.
- **Size** means U(N) body number: $X^{(k)}$ is the component of $X\in\mathrm{End}(\Lambda^F)$ in the irrep $V_k$
  (adjoint Casimir $2k(N+1-k)$). Then $f_k(H)=\|H^{(k)}\|^2/\|H\|^2$ and $w_k=\|P^{(k)}\|^2/d$ for the BPS projector.
- $P_n$ is the BLY chord-number distribution of the SUSY Hartle–Hawking state at $\lambda=2p^2/N'$, $j=0$.

### Step 1: $H$ is at most $(p-1)$-body [derived]
Write $I=A\sqcup S$ and $J=B\sqcup S$ with $S=I\cap J$, $|S|=m$. Normal ordering gives, up to an overall sign,
$$
\psi_I\psi_J^\dagger=\pm\,\psi_A\psi_B^\dagger\prod_{s\in S}(1-n_s),\qquad
\psi_J^\dagger\psi_I=\pm\,\psi_B^\dagger\psi_A\prod_{s\in S}n_s ,
\qquad \psi_A\psi_B^\dagger=(-1)^{(p-m)}\psi_B^\dagger\psi_A .
$$
- **$m=0$.** $\psi_I\psi_J^\dagger=(-1)^{p^2}\psi_J^\dagger\psi_I=-\psi_J^\dagger\psi_I$ for odd $p$, so
  $\{\psi_I,\psi_J^\dagger\}=0$.
- **$m=1$** (checked by hand for $p=3$). With the signs in place the two terms combine as
  $\psi_A\psi_B^\dagger[(1-n_s)+n_s]=\psi_A\psi_B^\dagger$, a pure $(p-1)$-body operator.
- **General $m$.** The combination is $\psi_A\psi_B^\dagger\,[\prod_s(1-n_s)\pm\prod_s n_s]$, with the relative sign
  fixed so that the $m=0$ case vanishes. The top-degree term $\prod_s n_s$ cancels, so the bracket has degree
  $\le m-1$ in the $n_s$. Hence the body number is $\le(p-m)+(m-1)=p-1$.

The sign bookkeeping for general $m$ is sketched, not written out. The conclusion is FGMS's statement that $H$
contains terms with up to $2\hat q-2$ fermions (FGMS §2, after eq. 2.27), read in U(N) language, and it is verified
numerically:

| $p$ | $N'$ | $f_0$ | $f_1$ | $f_2$ | $f_3$ | $f_4$ | $f_{\ge5}$ |
|---|---|---|---|---|---|---|---|
| 3 | 10, 12 | 0.30, 0.38 | 0 | 0.70, 0.62 | 0 | 0 | 0 |
| 5 | 10 | 0.008 | 0 | 0.060 | 0 | 0.932 | 0 |
| 5 | 14 | 0.052 | 0 | 0.057 | 0 | 0.892 | 0 |

The zeros are exact, $\le10^{-28}$. At $p=5$ the 2-body part comes from $m=3$ overlaps.

### Step 2: size per chord unit [approx → conj]
Lin's derivation of $\bar n=\mathrm{size}/q$ (his eq. 58) inserts the size operator as a chord. Each $H$-chord crossing
it contributes the size of one Hamiltonian term. In the BLY super-chord theory one unit of chord number is one
$Q$–$Q^\dagger$ pair across the cut, i.e. one $H$. Transcribing Lin's argument therefore suggests
$$
\boxed{\;\text{the chord-number-}n\text{ sector of }P\text{ has U(N) size}\ \approx(p-1)\,n\;}
\qquad(\text{sizes truncated at }\min(F,N-F)).
$$
Consequences:
1. **For $p=3$ this is $k=2n$**, the empirical $w_{2n}\approx P_n$.
2. **The naive $k=p\,n$ is excluded.** For odd $p$ it would put the $n=1$ weight at odd $k$, which vanishes exactly
   at half filling (particle–hole; progress report §4.4).
3. **The one-unit sector should inherit $H$'s own size distribution.** It should not be concentrated exactly at
   $k=p-1$.

This is not a derivation for $\mathcal N=2$. Lin's argument is for Majorana DSSYK at $\lambda\to0$, where each chord's
size is sharp up to $O(\lambda)$ Poisson overlaps (his eq. 59). At our $\lambda=18/N'\approx1.3$–$2.3$, overlaps are
$O(1)$.

### Step 3: numerical test [num]
**Test 1: $p=5$, the case where $p-1\ne2$.** The prediction is weight at $k=4$ (plus $H$'s own $k=2$ admixture).

| $N'$ | $a$ | share of $1-a$ at $k=2$ / 4 / 6 | $H$: traceless share at $k=2$ / 4 | BLY $P_1/(1-P_0)$ |
|---|---|---|---|---|
| 10 | 0.9921 | 0.061 / 0.939 / – | 0.061 / 0.939 | 0.993 |
| 12 | 0.9740 | 0.067 / 0.926 / 0.007 | 0.068 / 0.932 | 0.984 |
| 14 | 0.9470 | 0.057 / 0.924 / 0.019 | 0.060 / 0.940 | 0.972 |

- The one-unit sector sits at $k=4$ and $k=2$ in exactly $H$'s proportions.
- The small $k=6$ weight tracks the $n\ge2$ chord share (0.016, 0.028), truncated because $k=8$ is not available.
- Two seeds agree to $\le10^{-3}$.

**Test 2: $p=3$, shapes conditional on nonzero size vs nonzero length.** The last bin absorbs truncated sectors.

| $N'$ | $w_{2n}/(1-a)$, $n=1,2,3$ | $P_n/(1-P_0)$, last = tail | $a$ vs $P_0$ |
|---|---|---|---|
| 8 | 0.931, 0.069 | 0.883, 0.117 | 0.771 vs 0.789 |
| 10 | 0.839, 0.161 | 0.807, 0.193 | 0.643 vs 0.671 |
| 12 | 0.735, 0.236, 0.028 | 0.727, 0.209, 0.064 | 0.526 vs 0.559 |
| 14 | 0.644, 0.283, 0.072 | 0.648, 0.247, 0.105 | 0.425 vs 0.459 |

- The $n=1$ conditional share converges to BLY's: the differences are $+0.048, +0.032, +0.008, -0.004$.
- So most of the 10–20 % termwise mismatch in the progress report's Table 3 is the **normalization** $a$ vs
  $a_\lambda=P_0$ (a finite-$p$ BPS-fraction effect), not the shape of the size distribution.

### Status
- **Derived:** $H$ is at most $(p-1)$-body for odd $p$ (sign bookkeeping for general $m$ sketched; verified
  numerically for $p=3,5$).
- **Numerical:**
  - the non-identity size weight of $P$ sits at multiples of $p-1$, with $H$'s own admixture ($p=3$, $N'\le14$;
    $p=5$, $N'\le14$);
  - conditional shapes match BLY's chord-sector shares.
- **[Superseded in part by D4, 2026-10-02: the bound $k\le(p-1)n$ is now exact at half filling, and $n\le1$ is exact for
  $p=3$.]**
- **Conjectural:** "size $\approx(p-1)\times$ chord number" as a statement about the super-chord Hilbert space. It is
  the $\mathcal N=2$, U(N)-covariant analogue of Lin's $\bar n=\mathrm{size}/q$. It needs an $\mathcal N=2$ version of
  Lin's eq. 58 in which $H$-chords, not single $Q$-chords, cross the size cut. Such a version would promote the
  table rows $w_{2n}\leftrightarrow P_n$ and $w_0=a\leftrightarrow P_0$ to a dictionary.
- **Related exact statement (progress report §4.3).** The decoder variance is a two-sided light-probe correlator in
  the BPS TFD, $T_2=\langle\langle P|U\otimes U^*|P\rangle\rangle/d$. Its U(N) orbit average measures mean size,
  $E[T_2]=1-4\langle k(N+1-k)\rangle_w/(N(N+1))$. This is the exact finite-$N$ analogue of Lin's eq. 53 with $s=1$.

---

## D4. Chord number = Wiener-chaos degree; exact selection rule $k\le(p-1)n$

*2026-10-02. Full note: `research/notes/chord_size_dictionary_2026-10-02.md` (derivations, tables, Figure 1, reproducibility).
Code: `src/chord_chaos.py`, `src/plethysm.py`, `scripts/run_chord_chaos.py`, `scripts/analyze_chord_chaos.py`,
`scripts/run_topsize_projection.py`, `scripts/run_chord1_sizes.py`; tests `tests/test_chord_chaos.py`.
Data: `results/data/chord_chaos_2026-10-01/`.*

**Definitions.**
- $P(\omega)$ is the BPS projector as a function of the Gaussian couplings. $P_{[n]}$ is its bidegree-$(n,n)$
  Wiener-chaos component, and $\pi_n=E\|P_{[n]}\|^2/d$.
- $J_{nk}=E\|P_{[n]}^{(k)}\|^2/d$ is the joint (chord, U(N) size) distribution, with marginals $\pi_n$ and $w_k$.

**Results.**
1. **Chord number = chaos degree [derived+verified].**
   - For OU-correlated couplings, $E\,\tau(P(\omega)P(\omega_x))=\sum_n\pi_nx^{2n}$ (Mehler).
   - In the chord expansion each contraction across the cut carries $x$, so this is BLY (4.5) with $q^\Delta\to x$.
2. **Selection rule [derived+verified].** $J_{nk}=0$ unless $k\le n(p-1)+1$; at half filling $k\le n(p-1)$.
   - Proof: invariance, then the $\mathrm{Sym}^n\Lambda^p$ plethysm for odd $p$ ($\max\ell(\mu)=n(p-1)+1$), then Koike's
     stable range, then particle–hole parity.
   - Verified: $J_{2,6}=0.0000(1)$ at $N'=12$; $J_{1,3}=0.0035$ at $N'=9$ (bound attained off half filling).
3. **Exact low sectors ($p=3$, half filling) [derived+verified].**
   - $J_{0k}=a\,\delta_{k0}$.
   - $P_{[1]}=\gamma H^{(2)}$ (Schur), so $J_{1k}=\pi_1\delta_{k2}$.
   - $\pi_1/\pi_0=m_1^2/(m_2-(1+M^{-1})m_1^2)$ $=7/27,\,24/55,\,220/351,\,52/63$ at $N'=8$–14. It tends to $2q/(1-q)$
     (BLY) as $p\to\infty$ at fixed $\lambda$, with corrections $\approx c(\lambda)/p^2$.
4. **Krylov form [derived].**
   - $\kappa_n=a\,p_n(0)^2$ (averaged DOS) exactly.
   - BLY's $P_n=a_\lambda p_n^{\rm DS}(0)^2$ at $j_R=0$.
   - $\kappa_n\approx\pi_n$ [numerical].
5. **Diagonal dominance [numerical].**
   - $J_{2,4}/\pi_2=0.47,0.71,0.82,\approx0.87$ ($N'=8$–14).
   - Top-size parts lie along $(H^n)^{(2n)}$, a lower bound saturated to 1–5%.
6. **Probe [derived].**
   - The chord-1 sector is weighted by $\hat\chi_2(U)$, not $x^2$ (the friend mode cancels in the probe's sign).
   - BLY's ≈2% variance agreement is a compensation of $+0.034$ ($n=0$) and $-0.028$ ($n=1$) at $N'=14$.

**Status.**
- **Exact:** items 1–4, 6.
- **Numerical:** item 5.
- **Open:**
  - closed form for $\pi_n$, $n\ge2$ (the finite-$N$ HH recursion);
  - an analytic smearing law;
  - whether $P_{[1]}\propto{:}H{:}$ holds exactly for $p\ge5$. It holds to 4 digits at $N'=10,12$.


## D5. The decoder fourth moment as a four-arc function: OU family, ℤ₂ chord rule, wormhole overlaps

*2026-10-05. Context: `research/notes/decoder_moments_2026-09-29.md` §11 (decomposition
$T_4=a^2+2a(T_2-a)+S+R$). Literature: BKS 1997, Pluma–Speicher 2019, BINT 2018, BBNR 2020, G. Lin 2026, CCMP 2026
(notes in `literature/notes/`). Normalization: $\tau=\mathrm{Tr}/d$ on $\Lambda^F$, so $T_{2m}(\mathbb 1)=1$.*

### Definitions
- **Couplings.** $C=(C_J)_{|J|=p}$ are i.i.d. standard complex Gaussians. $Q(C)=\sum_JC_J\Psi_J$, $P(C)$ is the BPS
  projector, $a=d/D$.
- **Correlated partner.** For a diagonal real correlation $K=(K_J)$, $|K_J|\le1$, let $\tilde C=KC+\sqrt{1-K^2}\,C'$ with
  $C'$ an independent copy.
- **Alternating functions.** $T_{2m}(K)=E\,\tau\big((P(C)P(\tilde C))^m\big)$.
- **Special cases.**
  - OU family: $K_J\equiv x$.
  - Decoder with probe $U_A=\prod_{i\in A}(-1)^{n_i}$: $K_J=\sigma_J=(-1)^{|A\cap J|}$, deterministic. This is exact
    because $U_AQ(C)U_A=Q(\sigma C)$, hence $U_AP(C)U_A=P(\sigma C)$.
  - Identity: $K\equiv1$. Independent models: $K\equiv0$.

### 1. One generating function for decoder and OU [derived, exact at finite N]
- Write $F(C)=P(C)\otimes P(C)$ and $\mathcal B(A\otimes B,\,C\otimes D)=\tau(ACBD)$, so that
  $T_4(K)=E\,\mathcal B(F(C),F(\tilde C))$.
- The multivariate Mehler formula for a diagonal correlation, valid also at $|K_J|=1$ by parity of Hermite
  polynomials, gives
  $$T_4(K)=\sum_\alpha K^\alpha\,t_\alpha,\qquad t_\alpha=\mathcal B(F_\alpha,F_\alpha).$$
  - $\alpha$ is a multi-index of chaos quanta in $(C,\bar C)$, and $F_\alpha$ is the corresponding Hermite component of
    $F$.
  - The same formula with $F\to P$ gives $T_2$.
- **Chord reading.** A quantum is a Wick contraction between an arc of flavour $C$ and an arc of flavour $\tilde C$,
  weighted $K_J$. Contractions between the two $C$-arcs (or the two $\tilde C$-arcs) are internal to $F_\alpha$ and
  carry weight 1.
- **Four-arc ℤ₂ rule.** On the circle with arcs $P,P',P,P'$:
  - chords between adjacent arcs carry $K_J$;
  - chords between opposite arcs carry 1.
- **Contrast with BINT.** BINT's random matter (crossed 4-pt function) instead weighs adjacent-arc chords $\tilde q$
  and opposite-arc chords $\tilde q^2$, a product rule. For a fixed ℤ₂ probe all three pairings of the four $U$'s give
  the same sign $(-1)^{\#U\text{ between the chord ends}}$, so there is no sum over $U$-pairings and BINT's R-matrix
  does not apply by substitution.
- **Consequences.**
  - The OU family and the decoder are **two evaluations of the same coefficients $t_\alpha$**: on the diagonal
    $K=x\mathbb 1$, and at the cube vertex $K=\sigma$.
  - Because the ensemble is $S_N$-invariant, the site average is
    $E_A[\sigma^\alpha]=\mathcal K_s(o(\alpha))/\binom Ns$, a Krawtchouk polynomial with $s=|A|$. Here $o(\alpha)$ is the
    number of sites covered an odd number of times by the index sets of the cross-arc chords, i.e. the support of the
    monomial $\prod\Psi$.
  - For a single site this is $1-2o(\alpha)/N$. For *independent uniformly random* index sets,
    $E[(-1)^{|A\cap\cup J|}]=x_s^{|\alpha|}$ exactly, with $x_s=\mathcal K_s(p)/\binom Ns$ (`x_exact` in the
    parity-family table).
  - **So the decoder equals OU at $x=x_s$, coefficient by coefficient, exactly when the cross-chord index sets are
    uncorrelated.**
  - **Double scaling [derived, leading order].** With $p^2/N$ fixed and a heavy probe $s\propto\sqrt N$, the overlap
    correlations shift $o$ by $O(1)$ and change the weights by $O(s/N)\to0$, so decoder ≡ OU($x_s$) for every
    $T_{2m}$. For a single site ($s=1$), $x_s\to1$ in strict double scaling and the decoder trivializes.
  - **Finite $p$ [num, D4 item 6].** Correlations are large. The chaos-1 sector $P_{[1]}=\gamma H^{(2)}$ has forced
    overlaps ("friends") and weight $\hat\chi_2(U)=0.505$ instead of $x_s^2=0.327$ at $N'=14$. This is the same effect
    as D3's size per chord unit $p-1$ rather than $p$.
  - **Hence the natural test (Phase A2):** compare the decoder with OU at the *matched* $x_*$ defined by
    $T_2^{\rm OU}(x_*)=T_2^{\rm dec}$, not at $x_s$.
  - **Leading-order prediction for $x_*$ [derived at chaos order ≤ 1; num].**
    - $T_2$ is dominated by chaos orders 0 and 1. The decoder weights chaos 1 by $\hat\chi_2(U)$ and OU by $x^2$.
      Hence $x_*\approx\sqrt{\hat\chi_2(U)}=\sqrt{1-8(N-1)/(N(N+1))}$.
    - Measured: 0.476 / 0.557 / 0.582 / 0.630 / 0.650 / 0.684 / 0.696 against predicted 0.471 / 0.537 / 0.588 /
      0.628 / 0.660 / 0.688 / 0.711 for $N'=8$–14. Within 2 %, except 3.7 % at $N'=9$ (odd, off half filling).
    - The naive $x_s$ is 0.25–0.54.
- **Remark.** At matched $T_2$, $T_4^{\rm dec}-T_4^{\rm OU}(x_*)=R^{\rm dec}-R^{\rm OU}(x_*)$, up to the $x$-dependence of
  $S^{\rm OU}$, which is a covariance of per-realization size data and is checked numerically. So the test needs no
  size decomposition.

### 2. Exact properties of the OU curve [derived]
- $T_2^{\rm OU}(x)=\sum_n\pi_nx^{2n}$ (D4).
- $T_4^{\rm OU}(x)=\sum_n\theta_nx^{2n}$ with $\theta_n=\sum_{|\alpha|=2n}t_\alpha$. Even powers only, since
  $P(-C)=P(C)$.
  - Unlike $\pi_n$, the $\theta_n$ need **not** be non-negative. $\mathcal B$ is not a positive form:
    $T_4=\tau(Z^2)$ with $Z=P(C)P(\tilde C)$, whereas the time-ordered $\tau(Z^\dagger Z)$ is a norm.
- $T_4(1)=1$ and $T_4(0)=c_0$, the two-model value.
- **Sum rule.** $\sum_n n\,\theta_n=2\sum_nn\,\pi_n$, since $\partial_xT_4=2\,\partial_xT_2$ at $x=1$ by $P^2=P$.
- $T_4-2T_2+1=\tau\big((PP_xP-P)^2\big)\ge0$, the second moment of $\sin^2$ of the principal angles. For the decoder
  this equals $4\tau\big((PP_oP)^2\big)$, with $P_o=(P-UPU)/2$ the part of $P$ odd under $U$.
- **Remainder.** Per realization $Q_4=S+R$ as in §11, with $Y_x=P_x-a$. Taking $E_{C'}$ at $x=0$ makes
  $E[|Y'\rangle\rangle\langle\langle Y'|]$ $U(N)$-invariant, which is orthogonal to $\delta\mathcal P$. So
  **$R^{\rm OU}(0)=0$ exactly**, and $R^{\rm OU}(1)=(1-a)^2-S$.
- **Small-angle expansion** $\epsilon=1-x\to0$, with real coupling coordinates $\omega_I$:
  - $PP_xP-P=-B^\dagger B+O(\delta\omega^3)$, where $B=\sum_I\delta\omega_I(1-P)\partial_IP\,P$.
  - Then $T_2=1-\frac{2\epsilon}{d}E\,\mathrm{Tr}\,g_{\rm tot}+O(\epsilon^2)$, so $\sum_nn\pi_n=E\,\mathrm{Tr}\,g_{\rm tot}/d$.
  - And
    $$T_4-2T_2+1=\frac{4\epsilon^2}{d}E\Big[\mathrm{Tr}\,g_{\rm tot}^2+2\sum_{I,J}\mathrm{Tr}\,g_{IJ}^2\Big]+O(\epsilon^3).$$
  - Here $G_{IJ}=P\partial_IP(1-P)\partial_JP\,P$ is the operator-valued quantum geometric tensor (CCMP A.7–A.11),
    $g_{IJ}=(G_{IJ}+G_{JI})/2$, and $g_{\rm tot}=\sum_Ig_{II}$.
  - The Berry curvature $F_{IJ}=G_{IJ}-G_{JI}$ **cancels** between the $IJIJ$ and $IJJI$ Wick orderings. Isotropic OU
    noise measures the variance of the quantum metric, not G. Lin's TOC − OTOC curvature moment.

### 3. Two-sided-wormhole reading [derived in double scaling (Majorana); N=2 by analogy]
- In the q-Gaussian GNS space, left and right multiplications commute. Hence
  $T_4^{\rm OU}(x)=\langle\Psi_{h_x,h}|\Psi_{h,h_x}\rangle$, with
  $|\Psi_{h,h'}\rangle=L(P(h))\,R(P(h'))\,\Omega$.
  - $|\Psi_{h,h'}\rangle$ is the zero-energy state of the two-sided chord Hilbert space whose left boundary evolves with
    $H(\omega)$ and right boundary with $H(\omega_x)$.
  - Chords are coloured by the boundary that emitted them. The inner product is the q-Fock product with
    $\langle h,h_x\rangle=x$.
- So **$T_4^{\rm OU}$ is the overlap of a mismatched BPS wormhole with its mirror**, and
  $T_2^{\rm OU}=\langle\Omega|\Psi_{h,h_x}\rangle$.
- At finite N this is just vectorization, $|P_1P_2\rangle\rangle=(P_1\otimes P_2^T)|1\rangle\rangle$. Its content is
  that the double-scaled limit is a two-colour chord computation (BBNR's oriented chords with an origin colour), not a
  matter-chord R-matrix.
- **Out of time order.** By BKS (§4, Def. 4.1), $(P(\omega),P(\omega_x))$ are values of the q-OU process at two times.
  $T_2$ is time-ordered and given by the q-Mehler kernel. $T_4$ has operator times $(s,t,s,t)$ and is not determined by
  the Markov kernel.

### 4. Decoder as a conjugation flow [derived]
- $U=e^{i\pi n_i}$. So $P(\theta)=e^{i\theta n_i}Pe^{-i\theta n_i}$ interpolates $P\to P'$, a CCMP conjugation
  deformation with $\Lambda=n_i$. Its couplings are $C_J\to e^{i\theta}C_J$ for $J\ni i$, which is measure-preserving.
- The quantum metric along it is $\hat\Lambda(1-\hat\Lambda)$, $\hat\Lambda=Pn_iP=P-O_T$. Hence
  $T_2(\theta)=1-(1-\cos\theta)\cdot2\tau(\hat\Lambda-\hat\Lambda^2)$ exactly. This is a two-term Fourier series,
  since $n_i$ has eigenvalues 0 and 1.

### 5. Free-liberation baseline [derived; Demni–Hamdi–Hmidi eq. 2]
- **Model.** Rotate $P$ against itself by a free unitary Brownian motion: $J_t=PY_tPY_t^*P$, $\tau(P)=a$, giving
  DHH parameters $\lambda=1$, $\theta=a$, $m_n(0)=1$. This is the Hilbert-space-isotropic null model at
  $D\to\infty$.
- **Moments.** With $u=e^{-t}$ and $b=1-a$, the moment ODE gives
  $T_2=a+bu$, so $\eta\equiv(T_2-a)/(1-a)=e^{-t}$ exactly, and
  $$T_4=a^2(2-a)+4ab^2u+(E-2ab^2t)u^2,\qquad E=1-a^2(2-a)-4ab^2 .$$
- **Remainder.** Relative to the free two-model value $c_0^{\rm free}=a^2(2-a)$ (i.e. $S_{\rm free}=a^2(1-a)$):
  $$R_{\rm lib}(\eta)=2a(1-a)(1-2a)\,\eta+\big[E+2a(1-a)^2\ln\eta\big]\eta^2 .$$
- **Checks.**
  - The ODE agrees with the closed form to $10^{-11}$ (`tests/test_ou_fourarc.py`).
  - $R_{\rm lib}(1)=1-2a+a^3$, which is $(1-a)^2-S_{\rm free}$, as it must be.
- **Sign structure.** The linear term changes sign at $a=\tfrac12$. $R(1)$ changes sign at the golden-ratio root
  $a=(\sqrt5-1)/2=0.618$ of $a^3-2a+1$. This is the threshold observed empirically in
  `decoder_moments_2026-09-29.md` §11.3.
- **Small angle.** $\kappa\equiv(T_4-2T_2+1)/(1-T_2)^2\to1/(1-a)$ as $t\to0$. This is the free-Poisson value for
  $B^\dagger B$, $B=(1-P)\,dW\,P$.
- **Non-analyticity.** The $\eta^2\ln\eta$ term is non-analytic. The OU curve is a power series in $x^2$, so OU can
  approximate liberation but not equal it.

### Status
- **Exact:** §1 (generating function, four-arc rule, Krawtchouk weights), §2, §4, §5.
- **Leading order in double scaling:** decoder ≡ OU($x_s$) for heavy probes; the wormhole-overlap reading (Majorana,
  via Pluma–Speicher and BKS; the N=2 version by analogy).
- **Tested (A1–A2, `research/notes/decoder_T4_fourarc_2026-10-05.md`) [num].**
  - The $x$-dependence of $S^{\rm OU}$ is $\le10^{-6}$.
  - The OU curve is close to free liberation (§5) at equal $T_2$. The deviation tracks the two-model non-freeness.
  - The decoder does **not** match OU at equal $T_2$. Its excess (localization) is $0.017$–$0.047$ in $T_4$ for
    $N'=8$–14, or 37 % → 19 % of $T_4^{\rm dec}-T_4^{\rm free}$.
  - The Gaussian site-local family reproduces it, as does its infinitesimal version (excess quantum-metric
    heterogeneity $\kappa$).

## D6. Exact uplift formula; decoder walls from the uplift map; the arithmetic origin of localization

*2026-10-12. Context: D5 and `research/notes/decoder_T4_fourarc_2026-10-05.md` (the decoder's fourth-moment excess over
isotropic deformations, "localization"). Scratch code for this section was promoted to `src/uplift.py`; see the research
note `research/notes/uplift_localization_2026-10-12.md` for numbers, commands and figures.*

### Setup
- Let e = the added mode (site i = N′−1) and "o" = the N = N′−1 old modes. Split
  Λ^F(N′) = Λ^F(o) ⊕ Λ^{F−1}(o)∧e and write Ψ = a + b∧e.
- The p-form splits as ω = ω₀ + η∧e. Here ω₀ is the old model's p-form (couplings not touching e) and η is a
  (p−1)-form on the old modes (couplings touching e).
- Q₀ = ω₀∧ and E = η∧. For p = 3, η is a 2-form, so E commutes with Q₀ (and E† with Q₀†).

### 1. Block form of the supercharge [derived]
- On sector F, Q(a + b∧e) = Q₀a + (Q₀b + εEa)∧e with ε = (−1)^F, using e∧a = (−1)^{deg a} a∧e.
- The BPS conditions QΨ = 0 = Q†Ψ read:
  - Q₀a = 0 and Q₀†b = 0;
  - Q₀b = −εEa;
  - Q₀†a = +εE†b.

### 2. Exact uplift formula [derived + verified]
- Hodge-decompose in the old model: a = h + Q₀x and b = g + Q₀†y, with h ∈ B_o^F and g ∈ B_o^{F−1} harmonic.
- Ea is Q₀-closed in a degree without cohomology (F+2 lies outside the BPS window), so it is exact. Likewise E†b is
  coexact.
- Solving the two inhomogeneous equations with Moore–Penrose pseudo-inverses gives a = h + εK′b and b = g − εKa, with
  K = Q₀⁺E (Λ^F → coexact Λ^{F−1}) and K′ = (Q₀†)⁺E† (Λ^{F−1} → exact Λ^F).
- Eliminating b:
  $$a=(1+K'K)^{-1}(h+\varepsilon K'g),\qquad b=g-\varepsilon K a,\qquad (h,g)\in B^0=B_o^F\oplus B_o^{F-1}.$$
- This is an explicit isomorphism from the uplift reference B⁰ onto the BPS space B of the N′-mode model. Because K′K
  maps harmonic forms into exact forms, (1+K′K)⁻¹ is lower-triangular in harmonic ⊕ exact. Hence:
  - h-states: a = h − (1+D)⁻¹K′Kh, with D = K′K restricted to exact forms;
  - g-states: a = ε(1+D)⁻¹K′g.
- **Verified:** relative residual outside the exact BPS space ≤ 4.4·10⁻¹⁵ and full rank d, at N′ = 8, 9, 10, 12, 13, 14,
  and for scaled site couplings s = 0.4, 2.5.
- **Relation to CV.** The CV decoder D_T, fidelity and dressing cost (CV §5) are defined abstractly. This formula
  makes the uplift explicit in terms of the old model's Green's operator Q₀⁺. To first order, the dressing an old BPS
  state h needs is b ≈ −εKh, so CV's first-order dressing cost is ‖Kh‖²/‖h‖².

### 3. First-order decoder and the walls [derived; num]
- **First order in η** (no mixing): h-states have decoder eigenvalue λ = 1/(1+σ²), and g-states λ = σ′²/(1+σ′²), where
  σ, σ′ are the singular values of K|_{B_o^F} and K′|_{B_o^{F−1}}.
- **Exact walls.** λ = 1 exactly iff η∧h = 0 (CV/Q3's protected states); λ = 0 iff ι_η g = 0.
- **[num] Near-wall states are first-order exact.**
  - At N′ = 14 the upper quantiles of the exact decoder (q ≥ 0.75) agree with 1/(1+σ²) to ≤ 3·10⁻³.
  - Decoder eigenvectors with λ ≥ 0.97 have ‖h‖² ≈ 0.99 and median ‖η∧h‖/‖h‖ ≈ 1.1, against 4–8 in the bulk.
  - Wall fractions (min(λ, 1−λ) < ε) from first order alone are 10–15 % below the exact ones.
- **[num] The bulk is not first-order:** first order underestimates T₂ by 5–12 % (N′ = 10–15; table in §5).
- **Dressing cost = time-integrated BPS two-point function [derived; tested].**
  - On Λ^{F+2}(o) the exact forms are invariant under H₀ = Q₀Q₀† + Q₀†Q₀, and there H₀ = Q₀Q₀†. Since Eh is exact,
    K†K = E†(Q₀Q₀†)⁺E gives
    $$\sigma^2\in\mathrm{spec}\,\mathcal D_h,\qquad \mathcal D_h=P\,E^\dagger H_0^{-1}E\,P=\int_0^\infty d\tau\;P\,E^\dagger e^{-\tau H_0}E\,P ,$$
    with P the projector on B_o^F. Likewise σ′² ∈ spec P′E H₀⁻¹E†P′ (P′ on B_o^{F−1}; on coexact forms H₀ = Q₀†Q₀).
  - Test: `tests/test_uplift.py::test_dressing_cost_is_inverse_hamiltonian_two_point_function` (agreement 10⁻⁹).
  - **[interp]** E = η∧ = Σ_{i<j} η_ij c_i†c_j† is pair creation with Gaussian couplings independent of the old
    model's: a random two-fermion matter operator (R-charge 2/q̂). So the first-order decoder spectrum is the
    spectrum (as an operator on the old BPS space, not just its trace) of the zero-temperature, Euclidean-time-
    integrated two-point function of this matter operator. The τ-integral converges because charge F+2 has no BPS
    states. This is the object a double-scaled matter-chord computation (BLY §5 type, with an H⁻¹ insertion) would
    have to produce.

### 4. Arithmetic origin of the walls [derived dimension count; num]
- Ea is exact in Λ^{F+2}(o). So E|_{B_o^F} maps a space of dimension b_F = dim B_o^F into the exact forms, of dimension
  e_{F+2} = C(N, F+2) − e_{F+5}, recursively.
- **Forced exact atoms per wall:** max(0, b_F − e_{F+2}). For the g-wall, use the coexact forms of Λ^{F−3}.
- **Near-wall states:** when the ratio r = e_{F+2}/b_F is close to 1, the restricted map is nearly square. Its small
  singular values pile up at a hard edge, which shows up as decoder eigenvalues pressed against the walls.
- For p = 3 the ratio is:

| N′ | 10 | 12 | 13 | 14 | 15 | 16 | 17 | 18 |
|---|---|---|---|---|---|---|---|---|
| r (h-wall) | 0.44 | 0.67 | 0.99 | 0.96 | 1.31 | 1.33 | 1.70 | 1.78 |

- So the hard edge sits at N′ = 13–14 and disappears from N′ = 15 on. Since b_F ∼ 3^{N/2} and e_{F+2} ∼ 2^N/√N, r grows
  exponentially. This is the spectral counterpart of CV §7.2's arithmetic threshold for exact uplift.
- **Exact atoms [num].**
  - The measured counts are 80, 160, 41, 82, 0 at N′ = 11–15, identical across realizations.
  - The refined bound gives 79, 158, 27, 54, 0. A small deterministic excess (+1, +2, +14, +28) remains unexplained.
- **Pre-registered prediction, confirmed [num].** From old-model (N′−1) data only, first order predicted the fraction
  of decoder eigenvalues within 0.001 of a wall to drop from 0.117 (N′ = 14) to 0.003–0.012 (N′ = 15). The exact N′ = 15
  values are 0.012–0.023 (6 realizations), against 0.128 at N′ = 14.

### 5. Three-operator closed form, orthogonal splitting, monotone self-dressing [derived + verified]

*Code: `OldModel.three_operators`, `decoder_blocks`, `decoder_closed_form` in `src/uplift.py`; tests in
`tests/test_uplift.py`; numbers from `scripts/run_uplift_blocks.py` → `results/data/uplift_2026-10-12/blocks.jsonl`.*

**Hodge pieces.** Write X = Λ^F(o) = H_X ⊕ Ex_X ⊕ Co_X and Y = Λ^{F−1}(o) = H_Y ⊕ Ex_Y ⊕ Co_Y (harmonic, exact,
coexact). By §2, a = h + x with x ∈ Ex_X, and b = g + z with z ∈ Co_Y. Define three operators:
- K_H = K|_{H_X} : H_X → Co_Y (first-order dressing of h);
- K′_H = K′|_{H_Y} : H_Y → Ex_X (first-order dressing of g);
- L = K|_{Ex_X} : Ex_X → Co_Y.

**Adjoint relation [derived].** Let x = Q₀w ∈ Ex_X and z = Q₀†v ∈ Co_Y.
- Using [E, Q₀] = 0: Lx = Q₀⁺EQ₀w = Q₀⁺Q₀Ew = P_Co(Ew). Hence ⟨z, Lx⟩ = ⟨Q₀†v, Ew⟩ = ⟨v, EQ₀w⟩ = ⟨E†v, x⟩.
- Using [E†, Q₀†] = 0: K′z = (Q₀†)⁺Q₀†E†v = P_Ex(E†v). Hence ⟨K′z, x⟩ = ⟨E†v, x⟩.
- So ⟨z, Lx⟩ = ⟨K′z, x⟩ for all such x and z, i.e.
$$P_{\rm Ex}\,K'\,P_{\rm Co}=L^\dagger .$$
Verified to 6·10⁻¹⁵ (N′ = 10–13).

**Exact uplift as an anti-Hermitian linear system.** The two equations x = εK′b and z = −εKa become
$$(1+\mathbb L)\begin{pmatrix}x\\z\end{pmatrix}=\begin{pmatrix}\varepsilon K'_Hg\\-\varepsilon K_Hh\end{pmatrix},\qquad
\mathbb L=\begin{pmatrix}0&-\varepsilon L^\dagger\\ \varepsilon L&0\end{pmatrix}=-\mathbb L^\dagger .$$
1+𝕃 is invertible for any coupling strength (spectrum 1 + iℝ). Moreover
(1+𝕃)⁻¹†(1+𝕃)⁻¹ = (1−𝕃²)⁻¹ = diag((1+R′)⁻¹, (1+R)⁻¹) with R = LL† on Co_Y and R′ = L†L on Ex_X.
- **[context; literature not checked]** The one-sided half, b = g − εKa (Q-closedness), has the form of the classical
  homological perturbation lemma for Q = Q₀ + (εE)·(e∧) with the Hodge homotopy Q₀⁺. The BPS (harmonic) condition
  adds the dual half, a = h + εK′b. We have not checked whether the two-sided anti-Hermitian form, or the orthogonal
  splitting below, appears in the literature.

**Orthogonal splitting [derived; verified ≤ 1.1·10⁻¹⁵, N′ = 9–15].** Let B_h = Ψ(B_o^F, 0) (the BPS states whose
n_e = 1 part is coexact) and B_g = Ψ(0, B_o^{F−1}) (those whose n_e = 0 part is exact). Then
$$B=B_h\oplus B_g,\qquad B_h\perp B_g .$$
- Proof: for Ψ_h, the source vector is (0, −εK_Hh); for Ψ_g it is (εK′_Hg, 0). Their images under (1+𝕃)⁻¹ have
  inner product ⟨u, (1−𝕃²)⁻¹v⟩ = 0, because (1−𝕃²)⁻¹ is block-diagonal. The harmonic parts h and g are orthogonal to
  everything else.
- Equivalent direct proof from the BPS equations alone: with b_h = Q₀†v and a_g = Q₀w,
  ⟨a_h, a_g⟩ = ⟨Q₀†a_h, w⟩ = ε⟨b_h, Ew⟩ and ⟨b_h, b_g⟩ = ⟨v, Q₀b_g⟩ = −ε⟨v, EQ₀w⟩ = −ε⟨b_h, Ew⟩. They cancel.

**Canonical uplift = minimum-norm uplift [derived; 2026-10-13].**
- In CV's cohomological sense, a BPS state uplifts the old class [h] if its empty-sector component lies in [h], i.e. its
  old-Hodge harmonic part is h. By §2 these states are exactly Ψ(h, g) for arbitrary g ∈ B_o^{F−1}; the freedom in g is
  CV's coset ambiguity (CV eq. (40)).
- By the splitting, ‖Ψ(h, g)‖² = ‖Ψ(h, 0)‖² + ‖Ψ(0, g)‖². So Ψ(h, 0) is the unique minimum-norm uplift, equivalently the
  one orthogonal to every uplift from the occupied sector.
- This resolves both ambiguities of CV §4.3 (harmonic representative and coset) using only the old model's metric.

**Closed-form decoder [derived; verified to 10⁻¹⁰].** In the (non-orthonormal) parametrization by (h, g) ∈ B⁰, the
decoder μ = PUP|_B (U = (−1)^{n_e}) is the generalized eigenproblem 𝒩v = μ𝒢v, with
$$\mathcal G=\begin{pmatrix}1+K_H^\dagger(1+R)^{-1}K_H&0\\0&1+K_H'^\dagger(1+R')^{-1}K'_H\end{pmatrix},$$
$$\mathcal N=\begin{pmatrix}1-K_H^\dagger(1-R)(1+R)^{-2}K_H&-2\varepsilon K_H^\dagger L(1+R')^{-2}K'_H\\
\text{h.c.}&-1+K_H'^\dagger(1-R')(1+R')^{-2}K'_H\end{pmatrix}.$$
- First order is L → 0: 𝒩, 𝒢 diagonal with σ², σ′².
- Relation to the integrated pair two-point functions 𝒟 = E†H₀⁻¹E on the closed forms of Λ^F(o), and
  𝒟′ = EH₀⁻¹E† on the coclosed forms of Λ^{F−1}(o):
  - K_H†K_H = P_H𝒟P_H;
  - K_H†L = P_H𝒟P_Ex;
  - R′ = P_Ex𝒟P_Ex;
  - K′_H†K′_H = P_H𝒟′P_H;
  - R = P_Co𝒟′P_Co.
- By push-through, A†f(ABA†)A = A†A f(BA†A). For example
  K_H†(1+R)⁻¹K_H = [𝒟(1 + P_Ex𝒟)⁻¹]_{HH} (checked to 3·10⁻¹⁵).
- So the **diagonal blocks are functionals of 𝒟 (h-block) and 𝒟′ (g-block) alone**. The h–g block additionally
  needs the operator K′|_{H_Y} itself, not only its Gram matrix: the link between the X- and Y-sides.
- First order keeps only the harmonic blocks P_H𝒟P_H and P_H𝒟′P_H.
- Physical reading of a dressing mode with singular value ℓ of L: the first-order cloud k = K_Hh (n_e = 1) is
  reduced to norm² ⟨k,(1+R)⁻¹k⟩. A fraction R(1+R)⁻¹ of it reappears as an exact component x of a (n_e = 0).

**Monotone self-dressing [derived].** Within B_h, 1 − λ is the generalized Rayleigh quotient
⟨k,(1+R)⁻²k⟩ / (‖h‖² + ⟨k,(1+R)⁻¹k⟩). This is ≤ ‖k‖²/(‖h‖² + ‖k‖²) for every h, because n ≤ m ≤ ‖k‖² and m ↦ m/(c+m)
is increasing. By Courant–Fischer, every ordered eigenvalue of the self-dressed h-block lies at or above its
first-order value 1/(1+σ²); likewise every g-block eigenvalue lies at or below σ′²/(1+σ′²). Verified (test). The
h–g mixing then adds 2‖C₂‖²_F/d ≥ 0 to T₂, where C₂ is the off-diagonal block in orthonormal bases.

**[num] First order → self-dressed (blocks exact, C₂ = 0) → exact (seeds 9000, 9001):**

| N′ | T₂ first → self → exact | mixing 2‖C₂‖²/d | H first → self → exact |
|---|---|---|---|
| 10 | 0.709 → 0.737 → 0.750; 0.689 → 0.731 → 0.747 | 0.013, 0.015 | 0.058 → 0.062 → 0.066; 0.057 → 0.063 → 0.067 |
| 11 | 0.637 → 0.663 → 0.671; 0.647 → 0.677 → 0.685 | 0.008, 0.008 | 0.0905 → 0.0905 → 0.0886; 0.0866 → 0.0867 → 0.0855 |
| 12 | 0.608 → 0.650 → 0.669; 0.633 → 0.674 → 0.687 | 0.019, 0.013 | 0.0832 → 0.0834 → 0.0878; 0.0867 → 0.0873 → 0.0911 |
| 13 | 0.566 → 0.606 → 0.617; 0.554 → 0.595 → 0.609 | 0.010, 0.014 | 0.1114 → 0.1104 → 0.1090; 0.1123 → 0.1122 → 0.1095 |
| 14 | 0.572 → 0.621 → 0.638; 0.575 → 0.627 → 0.642 | 0.016, 0.016 | 0.1047 → 0.1034 → 0.1075; 0.1038 → 0.1033 → 0.1067 |
| 15 | 0.525 → 0.572 → 0.587; 0.552 → 0.592 → 0.601 | 0.015, 0.009 | 0.1269 → 0.1247 → 0.1233; 0.1270 → 0.1241 → 0.1225 |

- Self-dressing raises T₂ by 0.026–0.052 and mixing adds 0.008–0.019. H moves by ≤ 0.003 under self-dressing
  (N′ ≥ 11) and by +0.004 (even N′) or −0.001 to −0.003 (odd N′) under mixing.
- A per-mode mean field reproduces the self-dressed T₂ to ≤ 0.005 (N′ = 11–14). It replaces R by
  ρ_j = ⟨k_j,Rk_j⟩/‖k_j‖² along each singular direction of K_H, giving
  μ_j = [1 − σ_j²(1−ρ_j)/(1+ρ_j)²]/[1 + σ_j²/(1+ρ_j)].
- Singular values of L: median 0.35–0.48, max 0.76–0.96 (N′ = 10–13), all < 1 at unit coupling ratio. This is
  **not** a theorem; ℓ scales with the size of η.
- Block spectra:
  - Even N′: spec M_g = −spec M_h per realization (to 7·10⁻¹⁵ at N′ = 10, 12, 14; the λ ↔ 1−λ symmetry acts
    block-wise), and μ_h ≥ −0.39.
  - Odd N′: μ_g ≤ −0.21, so every g-state sits on the n_e = 1 side.

### 6. The fourth moment: heterogeneity, not walls [identity; num]

*Script: `scripts/analyze_uplift_heterogeneity.py` → `results/data/uplift_2026-10-12/heterogeneity_analysis.json`,
`results/figures/uplift_2026-10-12/fig_heterogeneity.pdf`.*

**Identity [exact].** With y = μ² = cos²2θ,
$$T_4=T_2-H,\qquad H\equiv\langle y(1-y)\rangle=\tfrac14\langle\sin^2 4\theta\rangle .$$
So at equal T₂ the localization (D5) is L = T₄^dec − T₄^OU(x*) = H^OU(x*) − H^dec. **Localization is exactly the
decoder's deficit of weight at intermediate angles.**

**The wall hypothesis is refuted [pre-registered, num].**
- The hard edge collapses at N′ = 15 (§4). If localization were the wall pile-up, it would drop there.
- Measured: L = 0.0432 (N′=13), 0.0467 (N′=14), 0.0454 ± 0.0003 (N′=15; 8 samples, K = 8). The localization is a property
  of the bulk of the angle distribution.

**Prediction from the old model [num].** Take H^dec ≈ H_first, with
H_first = ⟨y(1−y)⟩ at y = ((1−σ²)/(1+σ²))² and the analogue for σ′. Then T₄ ≈ T₂ − H_first.

All-site decoder means from the OU campaign (D5 data; N′ = 15: `ou_fourarc_2026-10-05/N15_{a,b}`). H_first is the mean over 13 seeds (6 at N′ = 15). L = T₄^dec − T₄^OU(x*) (D5).

| N′ | T₄^dec | T₂ − H_first | error | liberation error | OU(x*) error | L measured | L predicted | H_first | H^dec |
|---|---|---|---|---|---|---|---|---|---|
| 8 | 0.7778 | 0.7884 | +0.0105 | -0.0151 | -0.0175 | 0.0175 | 0.0280 | 0.0291 | 0.0396 |
| 9 | 0.6666 | 0.6694 | +0.0029 | -0.0312 | -0.0293 | 0.0293 | 0.0322 | 0.0625 | 0.0653 |
| 10 | 0.6790 | 0.6869 | +0.0079 | -0.0315 | -0.0321 | 0.0321 | 0.0400 | 0.0587 | 0.0666 |
| 11 | 0.5799 | 0.5793 | -0.0006 | -0.0444 | -0.0380 | 0.0380 | 0.0374 | 0.0897 | 0.0891 |
| 12 | 0.5964 | 0.6023 | +0.0059 | -0.0446 | -0.0417 | 0.0417 | 0.0475 | 0.0838 | 0.0897 |
| 13 | 0.5101 | 0.5077 | -0.0025 | -0.0539 | -0.0432 | 0.0432 | 0.0408 | 0.1107 | 0.1083 |
| 14 | 0.5283 | 0.5323 | +0.0040 | -0.0541 | -0.0467 | 0.0467 | 0.0507 | 0.1045 | 0.1085 |
| 15 | 0.4524 | 0.4490 | -0.0034 | -0.0606 | -0.0454 | 0.0454 | 0.0420 | 0.1264 | 0.1230 |

- For N′ ≥ 9 the T₄ error is ≤ 0.008 (≤ 0.006 for N′ ≥ 11), against 0.03–0.06 for liberation alone and 0.03–0.047 for isotropic OU.
- N′ = 8, where 74 % of the decoder spectrum is exact atoms, is the outlier (+0.011).
- The localization is predicted to within 2–14 % for N′ ≥ 11 (25 % at N′ = 10).

**Why it works [num].**
- It is **not** an identity: first order underestimates T₂ by 5–12 %.
- Per-state overlap tracking (N′ = 12–14) maps first-order states onto exact eigenstates. y shifts up by about +0.10
  for y₁ < 0.6, falling to +0.005 at the walls, with a spread of sd 0.12–0.16. The per-state change of y(1−y) runs
  from +0.07 (y₁ < 0.1) to −0.045 (y₁ ≈ 0.65), and these changes cancel in the mean.
- The block decomposition (§5) splits the cancellation.
  - Self-dressing barely moves H. Its small effect drifts with size: +0.004 at N′ = 10, about 0 at N′ = 11–12, and
    −0.0025 at N′ = 15.
  - Mixing gives the parity-alternating residual (+0.004 at even N′, −0.002 at odd N′).
  - The drift is why the odd-N′ residual grows. H_first ≈ H^dec is an empirical approximation at these sizes, not an
    asymptotic statement.
- *Superseded:* an earlier reading, "first order misplaces states by y ↔ 1−y", is withdrawn. The folded distributions
  of min(y, 1−y) do agree (KS 0.02–0.03 for N′ ≥ 13, against 0.07–0.09 for y), but per-state tracking shows a shift
  and spread, not a reflection.

**Parameter-free variant [num].** Taking T₂ from the BLY chord formula (D2) instead of from the decoder:
errors in T₄ are +0.025, +0.013, +0.018, +0.003, +0.010, −0.006, +0.002, −0.012 at N′ = 8–15 (≤ 0.012, i.e. ≤ 2.7 % of T₄, for N′ ≥ 11). These errors are dominated by the 1–2 % error of the chord T₂ (D2). No input from the N′-mode model is used.

**Out-of-sample test at N′ = 16 [pre-registered, `results/data/uplift_2026-10-12/prereg_N16.json`].**
- **Pre-registration.** Written 2026-10-12 21:58 from the N = 15 old model only (first-order seeds 9000–9002) and the
  N′ = 10–15 calibration, before any exact N′ = 16 computation. Exact runs: seeds 9000 and 9001, peak 10.8 GB each
  (memwatch cap 13 GB), 17–22 min each.
- **P1 (the main prediction): H_exact − H_first ∈ [0, +0.005] per seed.** Result: +0.0030 and +0.0042. **Pass.**
  Equivalently, T₄ = T₂ − H_first is within +0.003 / +0.004 of the exact T₄ (0.4829, 0.4887). This continues the
  even-N′ sequence +0.0079, +0.0045, +0.0034 (N′ = 10, 12, 14).
- **P2 (no exact atoms):** pass.
- **P5 (T₂ ∈ 0.585–0.615):** 0.6075 and 0.6119; pass.
- **P3 (wall fractions 5–20 % above first order):** ε = 0.1 and 0.001 pass. ε = 0.03 and 0.01 *narrowly fail*: the
  exact/first ratios are 1.19–1.26, so the walls are slightly more enhanced over first order than calibrated.
- **P4 (folded distributions):**
  - KS(min(y, 1−y)) = 0.018 and 0.025 (≤ 0.035) and KS(y) = 0.088 and 0.107 (≥ 0.05); pass.
  - The fraction of y < ½ is 0.367 (pass) and 0.358 (*narrow fail* of the 0.36–0.40 window).

### Status
- **Exact:**
  - §§1–2;
  - the atom characterization and the dressing-cost identity (§3);
  - the adjoint relation, the anti-Hermitian linear system, the orthogonal splitting B = B_h ⊕ B_g, the closed-form
    block decoder, and monotone self-dressing (§5);
  - T₄ = T₂ − H (§6).
- **Numerical:**
  - first-order accuracy near the walls;
  - the hard-edge collapse at N′ = 15 (pre-registered);
  - H_first ≈ H^dec to 1–7 % for N′ = 11–15 (12 % at N′ = 10), which predicts T₄ to ≤ 0.008 and the localization to 2–25 %
    (liberation alone misses T₄ by 0.03–0.06);
  - the pre-registered N′ = 16 test: the main prediction (H residual, i.e. T₄ to +0.003/+0.004) passed; two side
    predictions (wall fractions, y < ½ fraction) narrowly missed.
- **Refuted:** "localization = wall pile-up" (pre-registered drop at N′ = 15 not seen).
- **Open:**
  - why first order preserves H (an approximate cancellation; no symmetry found);
  - whether the odd-N′ residual (+0.002, −0.001, −0.003, −0.004 at N′ = 9, 11, 13, 15) keeps growing;
  - the residual atom excess (§4);
  - the chord (BLY matter-chord with H⁻¹ insertion) computation of the spectra of the harmonic blocks of 𝒟, 𝒟′.
    That would make the first-order spectrum, hence H_first and the T₄ prediction, analytic.

---

## D7. The first-order dressing spectrum from chords: closed-form BPS→continuum two-point function, free compound Poisson law, and the BPS heterogeneity

*2026-10-14. Context: D6 §3 (dressing cost σ² ∈ spec 𝒟, 𝒟 = P E† H₀⁻¹ E P) and the open item "chord computation of the
spectra of 𝒟, 𝒟′". Code: `src/dssyk_n2.py` (closed forms), `src/chord_transfer.py` (transfer-matrix engine, independent
check), `src/chord_dressing.py` (finite-N machinery); tests `tests/test_chord_dressing.py`; data and logs
`results/data/chord_dressing_2026-10-14/`. Report: `research/tex/progress_2026-10-14.tex`.*

### Setup and conventions
- One-sided N=2 chord Hilbert space at fixed charge s (BBNR §§4.4–5.1): bosonic states |∅⟩, |n,X⟩ = XOXO…XO,
  |n,O⟩ = OXOX…OX; fermionic |n+½,X⟩ = XOX…X, |n+½,O⟩ = OXO…O. Q opens O with q^s and closes X; Q† opens X with q^{−s}
  and closes O. Inner product BBNR (4.37). λ = 2p²/N, q = e^{−λ}, chord charge s = (degree − N/2)/p.
- Matter (BBNR §7.1): a chiral single-chord operator of p′ fermions shifts the charge of the region it encloses by
  δ = p′/p = 2Δ and inserts q̃^N̂ at the junction, q̃ = q^δ. For the pair operator E (p′ = 2, p = 3): Δ = 1/3, δ = 2/3.
- Dressing problems (old model N = N′−1, F = ⌊N′/2⌋):
  - h family: s_out = (F − N/2)/p, s_in = s_out + δ (exact forms of degree F+2);
  - g family: the mirror s → −s (coexact forms of degree F−3).
  - Even N′: (s_out, s_in) = (1/6, 5/6) for both. Odd N′: h (0, 2/3), g (1/3, 1).
- Calibration (no fit): energy scale κ from the exact first moment of H₀ in the inner sector,
  κ = 2[C(d,p) + C(N−d,p)]/(q^{s_in} + q^{−s_in}); τ = 0 value c₀ = 2C(N−F,2) (h), 2C(F−1,2) (g), exact.

### 1. BPS→continuum two-point function in closed form [derived + verified]
- **Ground state.** QΨ = Q†Ψ = 0 for Ψ = |∅⟩ + Σ_m (x_m|m,X⟩ + o_m|m,O⟩) gives, for m ≥ 0,
  q^s o_m + q^{−1}o_{m+1} − q^m x_{m+1} = 0 and q^{−s}x_m + q^{−1}x_{m+1} − q^m o_{m+1} = 0, solved by
  x_m = q^{3m/2} H_m(−cosh λ(s+½) | q²)/(q²;q²)_m (BLY 3.6–3.7).
- **Matter overlap.** For the B-family basis b_m = Q|m+½,X⟩ of the inner sector, write Q†_in = q^{−δ}(open part of Q†_out)
  + (close part) and use Q†_out Ψ = 0:
  $$\langle b_m|\tilde q^{\hat N}|\Psi\rangle=(1-\tilde q q^{\delta})\,(\tilde q/q)^m\,(q^2;q^2)_m\,x_m .$$
- **Eigenvectors.** BBNR (5.5) gives the B-family recursion. With the eigenvalue convention of BBNR (5.3),
  Λ_s(φ) = 2q^{−1/2}[cosh λ(s−½) − cos φ], the eigenvector coefficients are α_m(φ) = q^{m/2}H_m(−cos φ|q²)/(q²;q²)_m.
  BBNR (5.10) writes H_m(+cos φ), which belongs to +cos φ in Λ. This is harmless for BBNR (symmetric measure), but it
  reflects φ → π − φ in any kernel with one BPS end. Checked directly against (5.5).
- **q-Mehler** then gives A(φ) := ⟨v_in(φ)|q̃^N̂|Ψ⟩/⟨∅|Ψ⟩ = (q^{2δ};q²)_∞ / (q^{s_in+3/2}e^{±iφ}, q^{½+δ−s_out}e^{±iφ}; q²)_∞.
- Combined with the B-family density ρ_B = (q², e^{±2iφ};q²)_∞ q^{−s_in}/(2πΛ) (BBNR 5.11), and using
  1/Λ = q^{s_in}/|1 − a e^{iφ}|² with a = q^{s_in−½} together with q̃q^{3/2+s_out} = a q², one factor cancels:
  $$G(\tau)=\int_0^\pi d\phi\, f(\phi)\,e^{-\tau\kappa\Lambda_{s_{in}}(\phi)},\qquad
    f=\frac{(ab;q^2)_\infty\,(q^2,e^{\pm 2i\phi};q^2)_\infty}{2\pi\,(ae^{\pm i\phi},be^{\pm i\phi};q^2)_\infty},\quad
    a=q^{s_{in}-1/2},\ b=q^{1/2+\delta-s_{out}} .$$
  **f is the normalized Al-Salam–Chihara measure** (base q², x = cos φ).
- **Checks.**
  - Transfer-matrix engine: moments Λ^k for k = −2…3 agree to 1e−15 at several (λ, s, δ).
  - For s_in < ½ the measure acquires a point mass at zero energy, (q^{1−2s_in};q²)/(q^{1−2s_out};q²). This equals the
    BPS→BPS weight of BLY (4.14) (to 1e−12).
  - Triple scaling (φ = λk, λ → 0): f ∝ k sinh πk |Γ(Δ − ¼ + s_out/2 + ik/2) Γ(Δ + ¼ − s_out/2 + ik/2)|²,
    which is LMRS (43) with j′ = s_out, j = s_in.
- Literature status: BLY give BPS→BPS; BBNR (7.17) give charge-integrated ground-state/continuum moments. We have not
  found the fixed-charge measure, the ASC identification or the time integral in BBNR, BLY or LMRS.

### 2. Time-integrated two-point function [derived + verified]
- At x_a = (a + a⁻¹)/2 the ASC polynomials collapse to p_n = (ab;q²)_n/(2a)^n. The Casoratian sum of the Stieltjes
  transform then gives
  $$\int_0^\infty G(\tau)\,d\tau=\frac1\kappa\int\frac{f}{\Lambda}=\frac1\kappa\,\frac{q^{s_{in}}}{1-q^{2\delta}}\;
    {}_2\phi_1\!\left(q^2,q^2;q^{2+2\delta};q^2,q^{2s_{in}-1}\right).$$
  Verified against quadrature to 1e−15.
- Mean first-order dressing cost: t₁ = c₀ × (this). As λ → 0,
  ∫f/Λ ≈ π/(2 sin πδ) λ^{δ−2}(2s_in − 1)^{δ−1}: the Schwarzian region dominates.

### 3. Higher moments: matter-chord diagrams [derived from the chord rules for uncrossed diagrams; crossed ones measured]
- E_η tr 𝒟^m is a sum over the m! Wick pairings (E_i ↔ E†_{π(i)}). Each pairing is a disk with alternating zero-temperature
  segments (charge s_out) and ∫dτ segments (charge s_in).
- **Uncrossed pairings** (π ∈ NC(m), Catalan many):
  - each face carries a single energy (BINT/BLY uncrossed rules);
  - faces with zero-temperature segments sit on the unique zero mode of s_out (rank one ⇒ "factorization");
  - a G-face with k segments gives (1/b) Tr S^k with S diagonal in energy, eigenvalue g(φ) = (b c₀/D_in) A(φ)/(κΛ(φ)).
  - Hence the free cumulants are κ_k(𝒟) = (1/b) Tr S^k: **𝒟 is free compound Poisson** with rate r = e/b and jump law g
    under the normalized B-family density. Atoms 1 − r at r < 1 reproduce the rank bound (D6 §4).
- **Crossed pairings** need the N=2 matter–matter R-matrix (not available; BLY §5 set up but did not evaluate).
  Measured instead with an unbiased three-vector estimator: the crossing term of E tr 𝒟³ is −0.44 %, −0.16 %, ≈ +0.1 % of
  the third moment at N′ = 10, 12, 14 [num].
- Exact finite-N identities [derived]:
  - E_η tr 𝒟 = (1/b) Tr S = (1/b) Tr 𝒟̄;
  - E_η tr 𝒟² = tr 𝒟̄² + (1/b) Tr S²;
  - the separable (Kronecker) model X = S^{1/2} Z C^{1/2} with C = 𝒟̄/t₁ reproduces both.

### 4. Finite-N tests and error budget [num] (p = 3, N′ = 10–16; 2–4 realizations; no fitted parameters)
- S is diagonal in the inner energy basis to 98.5–99.8 %.
- t₁ (h): the chord value is low by 4–15 %. The exact charged spectrum is broader than the double-scaled one
  (CV² 0.36–0.41 vs 0.32–0.36) and extends below the chord gap. No λ calibration removes this.
- W^(1): the chord prediction exceeds the exact first-order value by 4–9 %. Error budget:
  - planar law with exact inputs: +5.5 % → +3.5 % from N′ = 10 to 16. The exact spectrum has a softer lower edge; at
    N′ = 10, giving each inner level its own η restores the free edge, so it is a finite-number-of-couplings effect;
  - chord inputs: +1–3 %;
  - BPS factorization: < 1 %.
- T₂^(1): the chord value is low by 2–6 %.

### 5. BPS heterogeneity of the averaged dressing cost [derived identity; num]
- 𝒟̄ = E_η 𝒟 is not ∝ 1 at finite N.
  - The relative spread of its eigenvalues is 0.22–0.42 (h), 0.12–0.16 (g), roughly flat in N′ = 10–16.
  - K(τ) = P E_η[E† e^{−τH₀} E] P has zero spread at τ = 0 (Σ_I ψ_I ψ_I† = C(N−F,2) on degree F), and the spread grows
    with τ: about 1.5 %, 4 %, 12 %, 30 % and 50–75 % at τE_med = 0.1, 0.3, 1, 3, 10.
  - The lowest 5 % of inner levels carry 31–89 % of the variance.
- **Sum rule** [derived + verified to 1e−14]. Let c_i = ι_i and E_i = {Q, c_i} = (ι_i ω)∧ (pair creation with the old
  couplings C_{i··}). For harmonic α, Q† c_i α = 0, so c_i α = harmonic + coexact, and the coexact part has squared norm
  ⟨α|E_i† G E_i|α⟩. Summing over sites,
  $$B_c+\sum_i P E_i^\dagger G E_i P = F\,P,\qquad B_c=\sum_i P\,c_i^\dagger P_{F-1}\,c_i\,P .$$
  𝒟̄ is 72–92 % (traceless-)correlated with Σ_i P E_i† G E_i P, hence anti-correlated with the single-fermion BPS→BPS
  weight B_c. The heterogeneity the chords miss is therefore the microstate dependence of BPS-projected single-fermion
  two-point functions.
- **Incoherent and coherent parts** [num]. Write 𝒟̄ = Σ_n M_n/E_n with M_n = 2Σ_I Pψ_I|n⟩⟨n|ψ_I†P (rank ≤ C(N,2)).
  - If the M_n were independently Haar-rotated, the eigenvalue variance would be Σ_n [tr M_n² − (tr M_n)²]/E_n² (free
    cumulants add). Call this the free, or incoherent, estimate.
  - Lowest 5 % of levels: observed/free spread = 1.00 (N′=10), 1.09 (11), 1.15 (12), 1.37 (13), 1.41 (14), 1.72 (15)
    for the h family. The coherent share of that variance, 1 − (free/obs)², is therefore 0, 0.17, 0.25, 0.48, 0.50, 0.66;
    for the g family it is 0.04, 0.27, 0.54 at N′ = 11, 13, 15.
  - The incoherent part shrinks as the number of contributing levels grows (N_eff = 21 → 365), while the total stays flat.
  - The whole is further reduced by the exact sum rule Σ_n M_n = c₀P, which makes the bins anti-correlated.
- [interpretation] The coherent part is a property of the BPS microstates (IR-heavy vs IR-light, i.e. small vs large
  single-fermion BPS weight). Whether it survives at large N is the same question as the large-N fate of localization.

### Status
- **Exact:** §1 (closed form; point mass = BLY 4.14; LMRS limit), §2, the finite-N identities of §3, the sum rule of §5.
- **Derived under the chord rules:** the free compound Poisson law (uncrossed diagrams + rank-one BPS faces).
- **Numerical:** §4 error budget; crossed diagrams ≤ 0.5 % of the third moment; the heterogeneity has an incoherent part
  (free estimate, shrinking with N′) and a coherent part (growing share, 0 → 0.66 at N′ = 10–15).
- **Pre-registered, not yet run:** the N′ = 17 test (`results/data/chord_dressing_2026-10-14/prereg_N17.json`, report §7).
- **Open:** the N=2 matter–matter R-matrix (crossed diagrams); the soft lower edge (shared couplings); the large-N fate of
  the coherent heterogeneity.

---

## D8. The two-flavor pair uplift: obstructed mapping cones, the protected tower, and the arithmetic of exact lifts

*2026-10-23. Context: the two-flavor analogue of D6 (Stage 1 of the monotone-vs-fortuitous study). Model: the
symmetrized two-flavor N=2 SYK of CCSY (arXiv:2412.06902, eq. 3.18), as studied by CV §6.2 (eqs. 82–90, Figs. 6–9).
Code: `src/two_flavor.py`, `src/two_flavor_uplift.py`; tests `tests/test_two_flavor.py`; data
`results/data/two_flavor_2026-10-23/` (`scripts/run_two_flavor_uplift.py`; analysis and figures
`scripts/analyze_two_flavor_uplift.py` → `summary.json`, `results/figures/two_flavor_2026-10-23/`). Realizations:
symmetrized seeds 1–5 (N = 2–6) and 1–2 (N = 7); generic seeds 1–3 (N = 2–6) and 1 (N = 7). Report:
`research/tex/progress_2026-10-23.tex`.*

### Setup and conventions
- Complex fermions ψ_i, χ_i (i = 1…N), written with explicit daggers. Supercharge (CCSY 3.18; CV 82):
  Q = Σ_{ijk} C_ijk ψ_i†ψ_j†χ_k with C symmetric in (j,k); here C = (G + G^{T_jk})/2, G i.i.d. complex Gaussian.
  Q has charge (+2,−1) under (N_ψ, N_χ); sectors are (a,b) = (N_ψ, N_χ).
- V = Σ_i ψ_i†χ_i† commutes with Q (CCSY 3.19). The tower V^n|Ω⟩ is BPS: Q-closed by [Q,V] = 0, Q†-closed through the
  su(2) of V, V†, N_V (CCSY 3.20–3.23). "Fortuitous" is the orthogonal complement of the tower in each sector (CV §6.2).
- Control: generic (HTZ) couplings, C i.i.d.; there V does not commute with Q and there is no tower.
- The uplift N → N+1 adds the pair (ψ_e, χ_f), both with coupling index N+1. Slots T = (t_ψ, t_χ), operator convention
  Ψ = Ψ₀₀ + ψ_e†Ψ₁₀ + χ_f†Ψ₀₁ + ψ_e†χ_f†Ψ₁₁, with Ψ_T in the old sector (A − t_ψ, B − t_χ) when Ψ is in the new sector (A,B).
- Old-model notation: B(a,b) is the BPS space, b(a,b) its dimension and P_H its projector. H₀ = {Q₀,Q₀†}, and
  G = H₀⁺ = ∫₀^∞ dτ e^{−τH₀} on the non-harmonic part. Exact and coexact forms are im Q₀ and im Q₀†.

### 1. Block form, identities and slot equations [derived + verified]
- **Block form.** Collecting the terms of Q_{N+1} with j = e, with k = f, and with both:
  $$Q_{N+1}=Q_0+\psi_e^\dagger E+K\chi_f+\psi_e^\dagger Z\chi_f,$$
  with
  - E = Σ_{i,k}(C_{e,i,k} − C_{i,e,k}) ψ_i†χ_k: even, charge (+1,−1);
  - K = Σ_{i<j}(C_{i,j,f} − C_{j,i,f}) ψ_i†ψ_j†: even, charge (+2,0);
  - Z = Σ_i (C_{e,i,f} − C_{i,e,f}) ψ_i†: odd, charge (+1,0).

  This is an operator identity, verified in both models (test).
- **Nilpotency.** Q_{N+1}² = 0, taken component by component in ψ_e† and χ_f, gives
  [Q₀,E] = 0 (ψ_e† part), [Q₀,K] = 0 (χ_f part) and **[E,K] = {Q₀,Z}** (ψ_e†χ_f part).
  - *Note added 2026-10-25 (D9 §1.1).* All four operators are polynomials in the mutually anticommuting set
    {ψ_i†, χ_k}. Hence [Q₀,E] = [Q₀,K] = [E,K] = 0 and {Q₀,Z} = 0 identically, and the last identity holds as 0 = 0
    (checked to 1e−15 in both models). Nothing below needs it to be non-trivial; in particular
    {Q₀,𝒵} = E P_H K − K P_H E (§9) is unaffected.
- **Protection identity.** Write V_{N+1} = V_N + ψ_e†χ_f†. The ψ_e† component of [Q_{N+1},V_{N+1}] = 0 gives
  $$[E,V_N]=K .$$
  The χ_f component gives [K,V_N] = 0, which holds automatically since both are built from creation operators.
  - Directly: [ψ_i†χ_k, V] = −ψ_i†ψ_k†, so [E,V] = −Σ_{i<k}(E_ik − E_ki)ψ_i†ψ_k†.
  - The identity is therefore K_ik = E_ki − E_ik, i.e. C_{N+1,i,k} = C_{N+1,k,i}: the (j,k) symmetry of the couplings
    at the new index.
  - Verified to 1e−16 in the symmetrized model. It fails for generic couplings (test).
- **Slot equations (S).** QΨ = Q†Ψ = 0, slot by slot. Define the slot currents
  J₁₀ := EΨ₀₀ − KΨ₁₁ + ZΨ₀₁ and J₀₁ := K†Ψ₀₀ + E†Ψ₁₁ + Z†Ψ₁₀. Then

  | slot | from QΨ = 0 | from Q†Ψ = 0 |
  |---|---|---|
  | 00 | Q₀Ψ₀₀ = −KΨ₀₁ | Q₀†Ψ₀₀ = −E†Ψ₁₀ |
  | 10 | Q₀Ψ₁₀ = J₁₀ | Q₀†Ψ₁₀ = 0 |
  | 01 | Q₀Ψ₀₁ = 0 | Q₀†Ψ₀₁ = J₀₁ |
  | 11 | Q₀Ψ₁₁ = −EΨ₀₁ | Q₀†Ψ₁₁ = K†Ψ₁₀ |

  Verified on every BPS state at N = 4, in both models (test).
- **Consequences.**
  - Ψ₀₁ is closed and Ψ₁₀ is coclosed.
  - The non-harmonic parts of Ψ₀₀ and Ψ₁₁ are fixed by Ψ₀₁ and Ψ₁₀:
    coexact(Ψ₀₀) = −Q₀⁺KΨ₀₁, exact(Ψ₀₀) = −(Q₀†)⁺E†Ψ₁₀, coexact(Ψ₁₁) = −Q₀⁺EΨ₀₁, exact(Ψ₁₁) = (Q₀†)⁺K†Ψ₁₀.
- **Filtration.** Q raises w := t_ψ − t_χ (01 → {00, 11} → 10) by 0 (Q₀), 1 (E, K) or 2 (Z).

### 2. Adding one mode with an obstruction [derived + verified]
**Theorem A.**
- *Setting.* Q₀ is nilpotent on a finite-dimensional graded space H_o, and a fermionic mode μ is added.
  - Creation type: Q = Q₀ + μ†E, with E even and [Q₀,E] = 0 (forced by Q² = 0).
  - Write Ψ = a + μ†b with a, b ∈ H_o.
- *(i) BPS conditions.* Ψ is BPS iff Q₀a = 0, Q₀†b = 0, Q₀b = Ea and Q₀†a = −E†b.
- *(ii) Existence and uniqueness.* Hodge-decompose a = h + x and b = g + z (h, g harmonic; x exact; z coexact).
  - A solution exists iff P_H E h = 0 and P_H E† g = 0.
  - It is then unique: z = Q₀⁺E(h + x) and x = −(Q₀†)⁺E†(g + z).
- *(iii) The BPS space.* Hence
  $$B\;\cong\;\ker\big(P_HEP_H\big|_{B_a}\big)\;\oplus\;\ker\big(P_HE^\dagger P_H\big|_{B_b}\big).$$
  - It splits orthogonally as B_h ⊕ B_g: no harmonic part in the μ-occupied slot, resp. in the μ-empty slot.
  - Ψ(h, 0) is the minimum-norm uplift of h.
- *Annihilation type*, Q = Q₀ + Kμ: the slot roles swap.
  - a = h + z (coexact) and b = g + x (exact).
  - The obstructions are P_H K g = 0 and P_H K† h = 0.

*Proof.*
- (i) Expand QΨ = Q₀a + μ†(Ea − Q₀b) and Q†Ψ = Q₀†a + E†b − μ†Q₀†b.
- (ii), existence. Q₀b = Ea needs Ea to be exact.
  - Ex = EQ₀w = Q₀Ew is exact. Eh is closed, since Q₀Eh = EQ₀h = 0, and closed = harmonic ⊕ exact.
  - So Ea is exact iff P_H Eh = 0. The coclosed equation is the mirror image.
- (ii), uniqueness. The pair (x, z) solves the anti-Hermitian system of D6 §5. Its proof uses only [Q₀,E] = 0 and the
  exactness of Ea, so it carries over verbatim, and so do the orthogonal splitting and the minimum-norm property. ∎

*Cohomological reading.*
- Q is the mapping cone of E: (H_o,Q₀) → (H_o,Q₀), whose long exact sequence is CV eq. (24) with connecting map δ = [E].
  The dimension formula in (iii) is **CV eq. (26)**.
- The new content is the harmonic (metric) refinement with δ ≠ 0: the explicit uplift map, the obstruction on harmonic
  sources, the orthogonal splitting and the canonical minimum-norm uplift.
- In the half-filling sector of the one-flavor model used in D6, δ = 0 because the target of [E] lies outside the
  fortuity window. In other one-flavor sectors δ can be nonzero; CV §4.2 show only that the incoming and outgoing
  connecting maps cannot both be nonzero.

*Verified* (`tests/test_two_flavor.py::test_theorem_A_explicit_step1`, N = 4, and scratch checks at N = 5):
- explicit solve at N = 4, 5 in both models: admissible sources give a BPS state with residual ≤ 1e−12, uniquely;
  obstructed sources give no solution (residual ≥ 1, i.e. O(1));
- the splitting holds to 1e−13.

### 3. The pair: two obstructed steps [derived + verified; the E₂ page is num]
- **The two steps.**
  - Step 1 adds ψ_e: creation type with E. The intermediate model M has χ_f frozen empty and Q_M = Q₀ + ψ_e†E.
  - Step 2 adds χ_f: annihilation type over M with K′ = K + ψ_e†Z. [Q_M, K′] = 0 follows from [Q₀,K] = 0 and
    [E,K] = {Q₀,Z}.
- **Dimension formula** (Theorem A twice):
  - dim B_M(A,B) = dim ker(P E P|_{B(A,B)}) + dim ker(P E† P|_{B(A−1,B)});
  - dim B_{N+1}(A,B) = dim ker(P K′† P|_{B_M(A,B)}) + dim ker(P K′ P|_{B_M(A,B−1)}).
  - Verified in every sector: 0 mismatches in 950 (symmetrized) and 570 (generic) sector records at N = 2–6, and in
    all 81 sectors at N = 7 (seed 1).
- **Four-slot injectivity [num].** A BPS state with vanishing old-harmonic part in all four slots is zero, so
  B_{N+1} embeds in ⊕_T B(slot T). Rank deficit 0 everywhere at N = 2–7 (both models). Every BPS state is reconstructed
  from its four harmonic slot parts to 1e−14.
  - A proof is open. The obstacle: the old-harmonic part of Ψ₀₀ is h₀₀ + P_H E†d, where d comes from the M-coexact
    dressing. So it differs from the step-1 source.
  - Injectivity holds, but the smallest bare fraction comes close to zero. It is 1.8e−6 (symmetrized, seed 1) and
    5e−5 (generic) in (3,3) and (5,5) at N = 7→8, against 1e−4 to 4e−3 at N = 5, 6.
  - *Update 2026-10-25 (§9.4):* this is a near-zero old gap making the dressing system nearly singular. Injectivity is
    generic but fails on a codimension-2 set of couplings, where pure-dressing BPS states appear.
- **The E₂ page is not the answer [num].** The slot complex built from old BPS states is
  H(01) → H(00) ⊕ H(11) → H(10), with d₁ = ([K], [E]) and (h₀₀, h₁₁) ↦ [E]h₀₀ − [K]h₁₁. Then:
  - d₁² = 0, because [E][K] − [K][E] = [{Q₀,Z}] = 0, and dim B_{N+1} ≤ dim E₂;
  - equality fails in exactly four sectors per N, identically in both models.
  - At N = 6 these are (2,3), (3,5), (4,2), (5,4), with E₂ = 92, 78, 78, 92 against 14, 0, 0, 14 BPS states. At
    N = 7 they are (2,3), (4,2), (4,6), (6,5), with E₂ = 91, 330, 330, 91 against 0, 239, 239, 0.
  - A differential d₂ ∼ [Z − EQ₀⁺K + KQ₀⁺E] pairs E₂^{w=−1}(A,B) with E₂^{w=+1}(A+2,B−1), always with the maximal
    rank (e.g. 9, 78, 91 at N = 5, 6, 7).
  - So the pair uplift is not computed by the slot complex alone. The two-step formula is exact and accounts for d₂.
  - *Update 2026-10-25 (§9.1):* the full transferred differential D = d₁ + d₂, with d₂ = P_H(Z − EQ₀⁺K + KQ₀⁺E)P_H,
    computes the cohomology from old-model data alone.
  - *(Corrects a remark made earlier in this project, that the spectral sequence degenerates at E₂: that holds for
    one mode, not for the pair.)*
- **Coupling independence of dimensions [num].** The BPS dimensions of the symmetrized and generic models coincide in
  every sector for N′ ≤ 8.
  - So the symmetrization creates no extra cohomology at these sizes: it makes one state per diagonal sector special.
  - This bears on CCSY footnote 18, which suggests that tower states "could easily be lifted" under non-symmetric
    perturbations. At these sizes, no BPS state is lifted.

### 4. The tower: protection by current cancellation [derived + verified]
- **Uplift of the tower.** From [E,V] = K, [K,V] = 0 and E|Ω⟩ = 0:
  $$E\,V^n|\Omega\rangle=n\,K\,V^{n-1}|\Omega\rangle .$$
  The new tower is V_{N+1}^n|Ω⟩ = V^n|Ω⟩ + n ψ_e†χ_f† V^{n−1}|Ω⟩ (CV 85). So its slots are Ψ₀₀ = V^n|Ω⟩,
  Ψ₁₁ = nV^{n−1}|Ω⟩ and Ψ₁₀ = Ψ₀₁ = 0, and both nonzero slots are old-harmonic.
- **Check of (S).**
  - (Q₁₀): J₁₀ = EV^n|Ω⟩ − nKV^{n−1}|Ω⟩ = 0. This is the protection identity.
  - (D₀₁): J₀₁ = K†V^n|Ω⟩ + nE†V^{n−1}|Ω⟩ = 0. This is equivalent to Q_{N+1}†V_{N+1}^n|Ω⟩ = 0 (CCSY 3.23), since
    Q_{N+1}†V_{N+1}^n|Ω⟩ = χ_f†J₀₁ once Q₀†V^m|Ω⟩ = 0 is used.
  - The other six equations hold trivially.
- **The tower is bare.** All its slot components are old-harmonic, so its Hodge dressing vanishes at every N and its
  bare fraction is exactly 1.
  - Slot weights: (N+1−n)/(N+1) in 00 and n/(N+1) in 11, from the norms (CV 86). Errors ≤ 4e−16.
  - Reading: the would-be dressing source J₁₀ of the 00 component (EΨ₀₀) is cancelled by the 11 component (KΨ₁₁). This
    current cancellation is the V Ward identity.
- **Two notions of dressing.** CV call this uplift "genuinely dressed" (κ₀₀ > 0), because relative to the 00 channel the
  state needs its 11 component. In Hodge terms it has no dressing. Both statements are correct. This document uses
  *bare* for the Hodge notion.
- **Not a decoder eigenvector.** The tower is not an eigenvector of CV's slot decoders P Π_T P: the residuals
  ‖PΠ_TP t − λt‖ are 0.12–0.43 for N = 5→6 through 7→8, and they are equal for T = 00 and 11.

### 5. Bare states: the bare map, its generic count, and the tower excess [derived count; num]
- **The bare map.** Bare BPS states are the kernel of the bare map β = [Q; Q†] restricted to ⊕_T B(slot T).
  With every Ψ_T = h_T harmonic, (S) reduces to:
  - Kh₀₁ = 0 in closed(A+2,B−1), and Eh₀₁ = 0 in closed(A+1,B−2);
  - E†h₁₀ = 0 in coclosed(A−2,B+1), and K†h₁₀ = 0 in coclosed(A−3,B);
  - J₁₀ = 0 in closed(A+1,B−1), and J₀₁ = 0 in coclosed(A−2,B).
  - J₁₀ is closed once the h₀₁ equations hold, because Q₀J₁₀ = [E,K]h₀₁ = EKh₀₁ − KEh₀₁. J₀₁ is the mirror case.
- **Generic count.**
  - Let n₀₁ = dim ker[K; E] on B(01) and n₁₀ = dim ker[E†; K†] on B(10). Generic position gives
    n₀₁ = max(0, b₀₁ − c₁ − c₂), and n₁₀ likewise.
  - Then bare_gen = max(0, b₀₀ + b₁₁ + n₀₁ + n₁₀ − c(J₁₀) − c(J₀₁)), with the c's the dimensions of the (co)closed
    targets listed above.
- **Result [num].** The law is
  bare = bare_gen + [one tower state, in diagonal sectors where bare_gen = 0]. It holds with 0 exceptions:
  - 400 sector records (symmetrized) and 240 (generic) at N = 2–6, provided the exact sub-kernel dimensions are used for
    n₀₁ and n₁₀;
  - the generic sub-kernel count fails by +1 only at N = 5: ker[E†;K†] on B(2,2) and ker[K;E] on B(2,1), together with
    their particle–hole mirrors ker[K;E] on B(3,3) and ker[E†;K†] on B(3,4). These failures occur in both models and
    are unexplained; they are absent at N = 4, 6 and 7;
  - the tower excess occurs at N = 6 in the central sectors (3,3) and (4,4), and at N = 7 in every diagonal sector
    (n,n), n = 1–7, where the generic count is 0;
  - in the generic model there is never an excess, and its arithmetic walls coincide with the symmetrized ones (N = 7:
    (1,0) 8, (2,0) 13, (2,1) 56, (3,1) 7, (3,2) 50 and mirrors, in both models).
- **The monotone/fortuitous statement at the level of bare states.** All bare fortuitous states are arithmetic walls,
  determined by dimension counting and the same in both models. The tower is the only bare state that the counting does
  not produce.

### 6. CV's channel decoders: constrained rooms, two counting laws, and the tower [derived rooms; num laws]
- **Rooms** (from (S); the analogue of CV §7.1). im D_T lies in a room R_T inside the old sector (a,b):
  - R₀₁ = ker Q₀;
  - R₁₀ = ker Q₀†;
  - R₀₀ = B(a,b) ⊕ Q₀⁺K(ker Q₀|_{(a,b−1)}) ⊕ (Q₀†)⁺E†(ker Q₀†|_{(a−1,b)});
  - R₁₁ = B(a,b) ⊕ Q₀⁺E(ker Q₀|_{(a+1,b)}) ⊕ (Q₀†)⁺K†(ker Q₀†|_{(a,b+1)}).
- **Two laws [num].**
  - (L1) rk D_T = min(d_new, R_T), where d_new = dim B_{N+1} of the target sector.
  - (L2) dim(B ∩ im D_T) = max(0, b + rk D_T − R_T).
  - Both hold with 0 failures in 1476 (sector, channel, seed) records (symmetrized) and 868 (generic) at N = 2–7.
- **Corollaries.**
  - (C1) U_T(fortuitous) = 1 ⟺ d_new ≥ R_T, including diagonal sectors, where the tower lies in im D₀₀ and im D₁₁.
  - (C2) U_T = 0 ⟺ (b − rk D_T) > 0 [protected channels 00 and 11 on the diagonal], resp. (f − rk D_T) > 0 [otherwise].
    Here f is the fortuitous dimension.
  - (C3) The value is intermediate otherwise.
- **Reproduction of CV Fig. 6 [num].**
  - The 0/1 pattern and the intermediate cells are reproduced for every column of CV's figure, N = 2→3 through 7→8
    (five seeds for N ≤ 6, two at N = 7). Intermediate values agree within realization scatter, e.g. (1,2)/00 at 3→4
    is 0.30 here vs 0.45 in CV, and (3,4)/00 at 7→8 is 0.17 vs 0.13.
  - The predicted category (1, (0,1) or 0) matches every symmetrized cell.
  - The generic model has one exception, at N = 5: (2,2)/01 and its mirror (3,3)/10. There B ∩ (im D)^⊥ has one more
    dimension than generic, in both models; this is another N = 5 coincidence, distinct from the sub-kernel one.
  - CV's qualitative pattern (00 fails for χ-heavy sectors, which lift through +ψ) is the inequality d_new ≥ R_T
    sector by sector.
  - Particle–hole symmetry pairs (a,b; 00) with (N−a, N−b; 11). Equal intermediate values in CV's own Fig. 6 (e.g.
    0.45/0.45, 0.65/0.65, 0.11/0.11) are this symmetry acting realization by realization, as here.
- **The tower in CV's channels.**
  - The tower lies in A₀₀ and A₁₁ for all N: it is not in generic position. Yet dim A_T still equals the generic count
    (the tower is one of its directions).
  - In the single-mode channels 10 and 01, the tower behaves like any other vector (CV Fig. 7).
- **Tower dressing costs [num law + derived corollary].**
  - κ_T(tower) attains CV's bound (eqs. 87, 88) exactly iff ker D_T = 0, i.e. iff the canonical completion (85) is the
    minimum-norm uplift: 0 failures in 255 (n, channel, seed) records at N = 2–7.
  - With (L1) this becomes κ₀₀ = n/(N+1−n) ⟺ d_new(n,n) ≤ R₀₀(n,n).
  - Together with (C1): **the tower's 00 cost locks onto CV's bound exactly when the 00 channel stops lifting every
    fortuitous state of the sector** (both at d_new = R₀₀ when equality holds).
  - This explains CV's observation that the cost "rises onto [the bound] only once N ≳ n+3": n = 1 from N = 4,
    n = 2 from N = 6, n = 3 from N = 7 (n = 4 still below at N = 7).
  - Sector (3,3) shows both halves: at N = 6→7 the cost is below the bound and U₀₀ = 1; at N = 7→8 the cost is on
    the bound and U₀₀ = 0.84.

### 7. First order: slot currents as time-integrated two-point functions [derived; num]
- **First-order dressing.** Replace Ψ_T by its harmonic source h_T on the right of (S). The six non-harmonic pieces are:
  - z₁₀ = Q₀⁺J₁₀(h) and x₀₁ = (Q₀†)⁺J₀₁(h);
  - z₀₀ = −Q₀⁺Kh₀₁ and x₀₀ = −(Q₀†)⁺E†h₁₀;
  - z₁₁ = −Q₀⁺Eh₀₁ and x₁₁ = (Q₀†)⁺K†h₁₀.

  κ⁽¹⁾ is the sum of their squared norms over ‖h‖². The first-order bare fraction is 1/(1+κ⁽¹⁾), as a generalized
  eigenproblem on the admissible sources.
- **Single-slot sources** (test: equal to the slot-diagonal-pseudo-inverse construction to 1e−8):
  $$\sigma^2_{00}=\langle E^\dagger G E\rangle+\langle K G K^\dagger\rangle,\quad
    \sigma^2_{11}=\langle K^\dagger G K\rangle+\langle E G E^\dagger\rangle,$$
  $$\sigma^2_{01}=\langle K^\dagger G K\rangle+\langle E^\dagger G E\rangle+\langle Z^\dagger G Z\rangle,\quad
    \sigma^2_{10}=\langle E G E^\dagger\rangle+\langle K G K^\dagger\rangle+\langle Z G Z^\dagger\rangle .$$
  - Each is a zero-temperature, Euclidean-time-integrated two-point function of a matter operator in the old model:
    flavor conversion ψ†χ (E), ψ-pair creation (K), or single ψ creation (Z).
  - G acts on (co)exact forms. For Z-terms, if J₁₀ is not closed the pseudo-inverse projects first.
- **Multi-slot sources interfere.** They share the currents J₁₀ and J₀₁, so κ⁽¹⁾ contains cross terms such as
  −2 Re⟨Eh₀₀, G Kh₁₁⟩ — mixed two-point functions of different matter operators.
- **The tower is the exactly dark state of this quadratic form.** J₁₀ = J₀₁ = 0 identically, while the incoherent
  (single-slot) sum is not small: 0.07–6.9 at N = 5 (seeds 1–5) and 0.07–1.8 at N = 6 (seeds 1–3).
  - The spread between realizations follows the old model's gap, which is as small as 0.03 at N = 5.
  - One exploratory realization with gap 0.0067 gave 60; it is not part of the production set.
- **[num] Accuracy.** First order reproduces the exact bare-fraction distribution. Sector means (sectors with ≥ 50
  states) agree to within 4.2 % at N = 5–6 and within 6.2 % at N = 7 (seed 1); first order is usually lower, e.g.
  N = 6 (3,3): 0.7014 vs 0.7322 and N = 7 (4,4): 0.6713 vs 0.7150. Walls coincide in every sector. That they must is
  clear for exact walls (κ⁽¹⁾ = 0); the converse is observed, not proved.
- **[num] Fortuitous states interfere too.** In the eigenbasis of the exact bare fraction, the median of
  (κ⁽¹⁾ − incoherent)/incoherent is −0.03 to −0.69 across sectors at N = 5–6 (seeds 1–5) and +0.00 to −0.81 at N = 7.
  - The ratio decreases with μ, from about +1 (constructive) for the most dressed states to −1 as μ → 1. At N = 7 a
    second, more destructive branch appears at small μ.
  - Near-complete cancellation is common among near-bare states.
  - **Correction of a working hypothesis:** "monotone = coherent, fortuitous = incoherent" is false. The tower is
    distinguished by *exact* cancellation at every N (an identity), not by the presence of interference.
  - Near-bare fortuitous states are tied to the hard edge of the bare map (§5), where bare_gen is close to
    threshold. The fraction with 0.99 < μ < 1 in the central sectors is 0.003 (N = 4), 0.07–0.08 (N = 5), 0.16 (N = 6:
    336 unknowns vs 339 constraints in (3,3), nearly square) and 0.065 (N = 7: 1092 vs 1428 in (4,4)). It peaks where
    the bare map is square. Partial N = 7 quantiles were seen before this was written, so it is a consistency check,
    not a blind test.

### 8. Monotone vs fortuitous: what the uplift distinguishes [summary]
| | tower (monotone) | fortuitous |
|---|---|---|
| Hodge dressing | 0 at every N (bare) | > 0, except arithmetic walls (peripheral sectors at N ≥ 7) |
| mechanism | J₁₀ = J₀₁ = 0 via [E,V] = K (V Ward identity) | partial interference |
| slot content | 00 and 11 only, weights (N+1−n)/(N+1), n/(N+1) | all admissible slots |
| bare count | the only excess over the generic count | equal to the generic count; same in the generic model |
| CV channels 00, 11 | always exact (non-generic position) | exact iff d_new ≥ R_T (generic position) |
| CV κ₀₀ | = n/(N+1−n) iff d_new ≤ R₀₀ | broad, growing (CV Fig. 8) |
| first-order dressing form | exactly dark | O(1), heavy-tailed |

### 9. Closed forms (added 2026-10-25) [derived + verified]

*Request: closed-form expressions for the pair uplift in old-model data (uplift map, obstruction and d₂, decoders,
injectivity). Code: `PairUplift.hpl_D`, `hpl_dim`, `cf_K`, `cf_bps`, `self_energies`, `reduced_system`, `cf_decoders`
in `src/two_flavor_uplift.py`; 4 new tests (20 in total); production task `closed` → `closed_N*_{sym,generic}.jsonl`;
counterexamples `scripts/find_pure_dressing_states.py` → `pure_dressing_*.json`.*

**Notation.**
- Slot vectors are Ψ = (Ψ₀₁, Ψ₀₀, Ψ₁₁, Ψ₁₀), ordered by w.
- Q_{N+1} = Q_d + Q_c:
  - Q_d = diag(−Q₀, Q₀, Q₀, −Q₀) is slot-diagonal (operator convention);
  - Q_c sends Ψ₀₁ ↦ KΨ₀₁ (00), EΨ₀₁ (11), ZΨ₀₁ (10); Ψ₀₀ ↦ EΨ₀₀ (10); Ψ₁₁ ↦ −KΨ₁₁ (10). It raises w by 1 (E, K) or 2 (Z).
- P_H is the slotwise old harmonic projector. E₁ = ⊕_T B(slot T) is the old BPS data. G = H₀⁺.

**9.1 Cohomology from old data: homological perturbation [derived + verified].**
- The data (ι, p = P_H, h = −Q_d⁺) form a strong deformation retract of (C, Q_d) onto (E₁, 0):
  Q_d h + h Q_d = ιp − 1, with h² = hι = ph = 0.
- Perturb by δ = Q_c. h preserves w and δ raises it, so (hδ)³ = 0, and the perturbation lemma (Brown–Gugenheim) holds
  with finite series.
- **Theorem B.** H(Q_{N+1}) ≅ H(E₁, D), where D = P_H Q_c (1 + Q_d⁺Q_c)⁻¹ P_H. Explicitly,
  $$D\eta=\big(0,\ P_HK\eta_{01},\ P_HE\eta_{01},\ P_H(E\eta_{00}-K\eta_{11}+\mathcal Z\eta_{01})\big),\qquad
    \mathcal Z=Z-E\,Q_0^+K+K\,Q_0^+E .$$
  - D = d₁ + d₂. Here d₁ is the slot complex of §3, and d₂ = P_H 𝒵 P_H : H(01) → H(10). There are no higher differentials
    (three filtration levels).
  - D² = 0. Its only component is P_H(E P_H K − K P_H E)P_H = P_H[E,K]P_H = P_H{Q₀,Z}P_H = 0, because (1 − P_H)K maps
    harmonic forms to exact ones and E keeps them exact.
  - Closed representatives (the HPL inclusion), for η ∈ ker D:
    Ψ_cl(η) = (η₀₁, η₀₀ − Q₀⁺Kη₀₁, η₁₁ − Q₀⁺Eη₀₁, η₁₀ + Q₀⁺(Eη₀₀ − Kη₁₁ + 𝒵η₀₁)).
  - **Dimension formula:** dim B_{N+1}(A,B) = dim E₁(A,B) − rk D(A,B) − rk D(A−2,B+1). It uses only old-model data and
    needs no intermediate model.
  - Verified: 0 mismatches in every sector record (950 symmetrized + 570 generic at N = 2–6, plus 81 + 81 at N = 7).
  - Reading: 𝒵 is a secondary, Massey-type operation. A bottom-slot class killed by d₁ reaches the top slot either
    directly (Z) or by zig-zags through the dressed middle slots (EQ₀⁺K, KQ₀⁺E). This sharpens §3: old cohomology plus
    [E], [K] is not enough, but old BPS data plus 𝒵 is exactly enough.

**9.2 BPS states in closed form [derived + verified].**
- **Theorem C.** Let 𝒦 = Q_d⁺Q_c + (Q_d†)⁺Q_c†. A slot vector Ψ is BPS iff
  (i) (1 + 𝒦)Ψ = P_HΨ, and (ii) P_H Q_cΨ = 0 = P_H Q_c†Ψ.
  The nonzero blocks of 𝒦 (row ← column):
  - 01 ← 00: −(Q₀†)⁺K†; 01 ← 11: −(Q₀†)⁺E†; 01 ← 10: −(Q₀†)⁺Z†;
  - 00 ← 01: Q₀⁺K; 00 ← 10: (Q₀†)⁺E†;
  - 11 ← 01: Q₀⁺E; 11 ← 10: −(Q₀†)⁺K†;
  - 10 ← 01: −Q₀⁺Z; 10 ← 00: −Q₀⁺E; 10 ← 11: Q₀⁺K.
- *Proof.*
  - (⇒) Decompose slotwise Ψ = h + r. Q_{N+1}Ψ = 0 gives Q_d r = −Q_cΨ. So Q_cΨ is exact, which gives P_H Q_cΨ = 0, and
    the coexact part of r is −Q_d⁺Q_cΨ. Dually, the exact part of r is −(Q_d†)⁺Q_c†Ψ. Hence r = −𝒦Ψ.
  - (⇐) From (i), Q_dΨ = −P_ex Q_cΨ, so Q_{N+1}Ψ = (P_H + P_co)u with u = Q_cΨ. Q_{N+1}² = 0 gives
    Q_d u = −Q_c(P_H + P_co)u. Induct on w: the lowest component of u is Q_d-closed, because Q_c raises w. Its coexact
    part therefore vanishes, and by (ii) so does its harmonic part. The next component is then closed, and so on, so
    Q_{N+1}Ψ = 0. Q_{N+1}†Ψ = 0 follows the same way by descending induction. ∎
- So if 1 + 𝒦 is invertible:
  - B_{N+1} = {(1 + 𝒦)⁻¹h : h ∈ 𝒜}, with 𝒜 = ker P_H(Q_c ⊕ Q_c†)(1 + 𝒦)⁻¹|_{E₁};
  - the four-slot harmonic map is injective.
- Verified (N = 2–6, both models):
  - every BPS state is reconstructed from its harmonic slot parts to ≤ 2.3e−13;
  - dim 𝒜 = dim B_{N+1};
  - the closed-form states have BPS residual ≤ 1e−13.

**9.3 Integrating out the middle slots [derived + verified].**
- **Middle slots.** Their dressing is explicit in terms of the outer slots: z₀₀ = −Q₀⁺KΨ₀₁, x₀₀ = −(Q₀†)⁺E†Ψ₁₀,
  z₁₁ = −Q₀⁺EΨ₀₁ and x₁₁ = (Q₀†)⁺K†Ψ₁₀.
- **Outer slots.** Substituting these leaves two "Dyson equations" for the outer components:
  $$P_{\rm Co}(1+\mathcal S_G)\Psi_{10}=Q_0^+\big(Eh_{00}-Kh_{11}+\mathcal Z\Psi_{01}\big),\qquad
    P_{\rm Ex}(1+\mathcal S'_G)\Psi_{01}=(Q_0^\dagger)^+\big(K^\dagger h_{00}+E^\dagger h_{11}+\mathcal Z'^\dagger\Psi_{10}\big).$$
  - 𝒮_G = E G E† + K G K† acts on slot 10, and 𝒮′_G = K† G K + E† G E on slot 01. Both are self-energies: Euclidean-time-
    integrated two-point functions of E and K through the old continuum.
  - 𝒵′ is 𝒵 one charge level down.
- **Reduced system.** For u = (z₁₀, x₀₁) ∈ Co(10) ⊕ Ex(01) this is 𝕄u = σ(h) with
  𝕄 = [[1 + 𝖠, −L], [−L′, 1 + 𝖠′]], 𝖠 = P_Co𝒮_G P_Co, 𝖠′ = P_Ex𝒮′_G P_Ex, L = Q₀⁺𝒵 and L′ = (Q₀†)⁺𝒵′†.
- **Structure.** 𝖠 and 𝖠′ are Hermitian PSD, and L′ = −L† + C with the leak
  $$C=\alpha_{00}\gamma_{00}^\dagger-\alpha_{11}\gamma_{11}^\dagger,\qquad
    \alpha_{00}=(Q_0^\dagger)^+K^\dagger|_{B_{00}},\ \gamma_{00}=Q_0^+E|_{B_{00}},\ \alpha_{11}=(Q_0^\dagger)^+E^\dagger|_{B_{11}},\
    \gamma_{11}=Q_0^+K|_{B_{11}} .$$
  B₀₀ and B₁₁ are the old BPS spaces of the middle slots.
- *Proof.*
  - From {Q₀, Z} = [E, K], Q₀Q₀⁺ = P_ex and Q₀⁺Q₀ = P_co one gets {Q₀, 𝒵} = Γ := E P_H K − K P_H E.
  - For x = Q₀w and z = Q₀†v this gives ⟨z, Lx⟩ = ⟨v, 𝒵Q₀w⟩ = −⟨L′z, x⟩ + ⟨v, Γw⟩.
  - Finally ⟨v, Γw⟩ = ⟨z, (γ₀₀α₀₀† − γ₁₁α₁₁†)x⟩. ∎
- **Reading of the leak.** γ₀₀h and α₀₀h are the first-order top-slot and bottom-slot dressings created by a
  middle-slot BPS state h. The leak is "top dressing → middle BPS state → bottom dressing". It vanishes when the middle
  slots carry no BPS states.
- Verified to 1e−13 (tests, and the scratch checks V1–V3).

**9.4 Uniqueness criterion: injectivity is generic, not universal [derived; num].**
- **Criterion.** The anti-Hermitian parts drop out of Herm 𝕄, leaving Herm 𝕄 = 1 + 𝖠 ⊕ 𝖠′ − ½[[0, C†],[C, 0]].
  - Hence Herm 𝕄 > 0 ⟺ ρ := ‖(1 + 𝖠′)^{−1/2} C (1 + 𝖠)^{−1/2}‖ < 2.
  - Then 𝕄 is accretive, hence invertible. 1 + 𝒦 is invertible too, since its middle rows are explicit, and the
    four-slot harmonic map is injective.
  - With no middle cohomology, C = 0, Herm 𝕄 ≥ 1 and ‖𝕄⁻¹‖ ≤ 1. This is the D6 contraction.
- **Numerics at unit coupling.**
  - N = 2–6: ρ < 2 in 426/450 (symmetrized) and 254/270 (generic) reduced systems. 𝕄 was invertible in all of them; the
    smallest σ_min(𝕄) is 0.31 and 0.11 respectively (N = 7: §9.8).
  - N = 7 (seed 1): σ_min(𝕄) drops to 0.013 in (3,3) and (5,5). The old model there has a near-zero mode (gap 0.0032 in
    sector (2,3)), so ‖𝒦‖ ≈ 148. This is the origin of the bare fraction 1.8e−6 noted in §3.
- **Counterexamples.** det 𝕄 is complex, so it vanishes on a codimension-2 set of couplings. Scaling two families of
  new-pair couplings and minimizing σ_min(𝕄) finds such points:
  - *generic model*, N = 4→5 (seed 2), sector (3,3): E-type couplings × 0.4216 and K-type × 12.58 give
    σ_min(𝕄) = 1.8e−16;
  - *symmetrized model*, N = 4→5 (seed 1), sector (2,2): C[N,·,·] × 7.280 and C[·,N,·] = C[·,·,N] × 12.30 give
    σ_min(𝕄) = 1.6e−16.
- **What happens there.** In both cases dim B_{N+1} is unchanged (65), yet one BPS state has bare fraction 0
  (injectivity deficit 1). It is a **pure-dressing** BPS state: all four of its slot components are non-BPS states of the
  old model. In the symmetrized example the tower of the same sector stays exactly bare.
- **Conclusion.** Four-slot injectivity (§3) is generic but false in general. The true statements are the criterion above
  and the always-valid two-step parametrization (§3).

**9.5 Admissibility in closed form [derived + verified].** For h ∈ E₁, Ψ = (1 + 𝒦)⁻¹h is BPS iff:
- (a) P_H K h₀₁ = P_H E h₀₁ = 0 and P_H E† h₁₀ = P_H K† h₁₀ = 0 (the w = 0 currents);
- (b) P_H[E h₀₀⁺ − K h₁₁⁺ + 𝒵 h₀₁] = 0, with h₀₀⁺ = h₀₀ + α₀₀†x₀₁ and h₁₁⁺ = h₁₁ + α₁₁†x₀₁;
- (c) P_H[K† h₀₀⁻ + E† h₁₁⁻ + 𝒵′† h₁₀] = 0, with h₀₀⁻ = h₀₀ + γ₀₀†z₁₀ and h₁₁⁻ = h₁₁ − γ₁₁†z₁₀.

Readings:
- (b) is Dη_Q = 0 for the HPL class representative η_Q = p_∞Ψ = (h₀₁, h₀₀⁺, h₁₁⁺, ·) of the BPS state; (c) is the dual
  statement.
- The same maps α, γ govern the leak and these corrections.
- Verified: (a)–(c) give the same subspace as Theorem C in every sector (N = 3–6, both models).

**9.6 Slot weights and decoder pencils [derived + verified].** For a BPS state with sources h and outer dressing u:
- ‖Ψ₀₁‖² = ‖h₀₁‖² + ‖x₀₁‖² and ‖Ψ₁₀‖² = ‖h₁₀‖² + ‖z₁₀‖²;
- ‖Ψ₀₀‖² = ‖h₀₀‖² + ⟨Ψ₀₁, K†GKΨ₀₁⟩ + ⟨Ψ₁₀, EGE†Ψ₁₀⟩;
- ‖Ψ₁₁‖² = ‖h₁₁‖² + ⟨Ψ₀₁, E†GEΨ₀₁⟩ + ⟨Ψ₁₀, KGK†Ψ₁₀⟩;
- κ := ‖Ψ‖²/‖h‖² − 1 = [‖u‖² + ⟨Ψ₁₀, 𝒮_GΨ₁₀⟩ + ⟨Ψ₀₁, 𝒮′_GΨ₀₁⟩]/‖h‖². This is the exact version of the first-order cost of
  §7, which drops the 𝕄 corrections and replaces Ψ₁₀, Ψ₀₁ by h₁₀, h₀₁.
- The bare-fraction decoder (μ) and CV's slot decoders (λ_T) are the generalized eigenvalues of (‖h‖², ‖Ψ‖²) and
  (‖Ψ_T‖², ‖Ψ‖²) on 𝒜. They agree with the exact spectra to 7e−13.
- **[num] Most of the exact correction is self-energy.** Three levels on the exact admissible sources (task `levels`):
  - first order: first-order sources, 𝕄 → 1;
  - self-dressed: exact sources, 𝕄 → 1 + 𝖠 ⊕ 𝖠′;
  - exact.

  Self-dressing moves the mean bare fraction 56–92 % of the way from first order to exact (N = 4–6, sectors with ≥ 50
  BPS states; median 0.67 symmetrized, 0.79 generic). Example: (3,3) at N = 6→7 gives 0.7014 → 0.7220 → 0.7322. The
  closed-form quadratic form reproduces ‖Ψ‖² to ≤ 6e−14.

**9.7 Special cases.**
- *One added mode.* There are no middle slots, so 𝖠 = 𝖠′ = 0 and C = 0, and 𝕄 = 1 + (anti-Hermitian). This is D6 §5.
- *The tower.* η = (0, V^n|Ω⟩, nV^{n−1}|Ω⟩, 0) has 𝒦η = 0 and σ(η) = 0, so u = 0 and Ψ = η, at every coupling (§4).

**9.8 N = 7 [num].** One seed of each model; 81 sector records per model.
- The dimension formula, the admissible dimension and the reconstruction all hold, with 0 mismatches and residuals
  ≤ 2e−11.
- The reduced system is invertible in all 49 sectors and accretive in 45 of them (per model). Its smallest σ_min is 0.013
  (symmetrized) and 0.051 (generic).
- Peak memory 5.5 GB.
- Totals at N = 2–7: the criterion ρ < 2 holds in 471/499 (symmetrized) and 299/319 (generic) reduced systems.

### Status
- **Exact:**
  - §1: the block form, the identities and the slot equations;
  - Theorem A, with its proof, splitting and minimum-norm uplift (§2);
  - the two-step dimension formula (§3);
  - the tower's protection, bareness and slot weights (§4);
  - the bare-map equations and the rooms (§§5–6);
  - the first-order formulas (§7);
  - *(added 2026-10-25, §9)* Theorem B: H(Q_{N+1}) ≅ H(E₁, d₁ + d₂), with d₂ = P_H(Z − EQ₀⁺K + KQ₀⁺E)P_H and a dimension
    formula in old-model data;
  - Theorem C: BPS states = (1 + 𝒦)⁻¹h with admissibility (a)–(c);
  - the reduced two-slot system, with self-energies EGE† + KGK†, K†GK + E†GE and the factorized leak;
  - the accretivity criterion ρ < 2 for injectivity;
  - the slot-weight formulas and decoder pencils.
- **Numerical laws** (N = 2–6 with 5 seeds, N = 7 with 1–2 seeds; generic control at all N):
  - four-slot injectivity at random couplings. It is *not* universal: it fails on a codimension-2 set of couplings, with
    explicit pure-dressing counterexamples (§9.4);
  - d₂ ≠ 0 in exactly four sectors per N, with maximal rank;
  - the room laws L1 and L2, and the tower-cost criterion;
  - bare = generic count + tower;
  - BPS dimensions are coupling-independent;
  - first order within 4–6 % for sector means.
- **Unexplained:** the one-dimensional coincidences at N = 5 in the near-central sectors (2,1), (2,2), (3,3), (3,4) — two
  particle–hole pairs — in both models.
- **Open:**
  - the geometry of the singular set of 𝕄 (pure-dressing states), and whether it approaches unit coupling as N grows;
  - simplifying the decoder pencils of §9.6 further, e.g. via push-through identities as in D6 §5;
  - the chord computation of the self-energies 𝒮_G, 𝒮′_G, the transfer 𝒵 and the leak maps α, γ (Stage 2). These are now
    the complete list of old-model inputs of the exact uplift.

---

## D9. Chords for the two-flavor uplift: the black-hole window, emission into it, and why the tower does not dress

*2026-10-25. Stage 2 of the monotone-vs-fortuitous study (Q12): a chord (double-scaled, "gravitational") description of the
two-flavor pair uplift of D8. Code: `src/two_flavor_chords.py` (dictionary), `scripts/run_two_flavor_chords.py` (tasks
window, blocks, resonance, pairing, tower, edgegap, twisted), `scripts/analyze_two_flavor_chords.py`; tests `tests/test_two_flavor_chords.py`;
data `results/data/two_flavor_chords_2026-10-25/` (`summary.json`); figures `results/figures/two_flavor_chords_2026-10-25/`.
Report: `research/tex/progress_2026-10-25.tex`. Realizations: symmetrized and generic seeds 1–3 (old N ≤ 6, N' ≤ 7) and
1–2 (N = 7, N' = 8–10), same couplings as D8.*

### What is imported and what is new
- **Imported.**
  - BBNR §§4.1–4.5: oriented chord diagrams, friend/enemy weights, and the localized fixed-charge transfer matrix.
  - BBNR (5.3), (5.25): continuum eigenvalues and the window |s| < 1/2 of the exact ground-state count.
  - BLY (3.17), (4.14): BPS fraction and BPS→BPS weight.
  - D7 §§1–2: the BPS→continuum spectral measure of a chiral matter chord (Al-Salam–Chihara) and its time integral.
  - CCSY §3: the model, the tower, the statement that symmetrization leaves the large-N G–Σ action unchanged (their §3.3),
    and the nonlinear charge constraint (their eq. 3.14).
- **Not in the papers we have.** A chord treatment of the two-flavor model. CCSY §3.3 and HTZ §3 use the G–Σ method only.
- **New here.**
  - The two-flavor chord dictionary (§1).
  - The window and its gapless edges in this model (§2).
  - The uplift as matter-chord emission, with its resonances (§3).
  - A chord-level account of the tower: twisted chords (§§1.2, 4).

### 1. The two-flavor chord dictionary [derived]
**1.1 Reduction to one-flavor chords.**

*Setting.*
- Let every term of Q create p_ψ ψ's and annihilate p_χ χ's: Q = Σ C_{IK} Ψ_I† X_K with |I| = p_ψ, |K| = p_χ. Our model has
  p_ψ = 2, p_χ = 1 and C_{ijk} − C_{jik} as the coupling of ψ_i†ψ_j†χ_k.
- Double scaling: p_ψ, p_χ → ∞ at fixed ratio, with λ := 2(p_ψ² + p_χ²)/N fixed, q = e^{−λ}.
- Sector (N_ψ, N_χ); A := N_ψ − N/2, B := N_χ − N/2.

*Claim.* The normalized moments m_k(A,B)/m_0(A,B) of H equal BBNR's fixed-charge chord moments at (q, s), with
  s = (p_ψA − p_χB)/(p_ψ² + p_χ²)   [= (2A − B)/5 for our model]. Q raises s by exactly 1.

*Derivation* (BBNR's steps, two flavors).
1. *Wick.* Moments are sums over oriented chord diagrams (BBNR 4.1–4.3). Each chord joins a Q to a Q† with the same
   (I, K).
2. *Pairs of chords.* For two chords the overlaps m_ψ = |I₁∩I₂| and m_χ = |K₁∩K₂| are independent, Poisson(p_ψ²/N) and
   Poisson(p_χ²/N).
   - *Enemy configurations* (BBNR 4.4 i, ii) put two ψ†'s (or two χ's) on a shared site with nothing between them. They
     therefore need m_ψ = m_χ = 0, with probability e^{−(p_ψ²+p_χ²)/N} = q^{1/2}.
   - *Friend configurations* gain a factor 2 per shared site of either flavor (BBNR 4.6–4.7). This gives
     E[2^{m_ψ+m_χ}] = q^{−1/2}.
   - *Crossings.* Ψ_I†X_K contains p_ψ + p_χ fermions (odd here), so crossings give (−1)^{#int}.
   - Hence m_k = 2^{−k} Σ (−1)^{#int} q^{(#e−#f)/2}, which is BBNR (4.10) with the combined λ.
3. *Chemical potentials.* Weight the trace with e^{−μ_ψ A_op − μ_χ B_op}.
   - A ψ-site of a chord carries ψ† at its O (Q) end; a χ-site carries χ there.
   - So at the base point of the trace, a right-pointing chord (O before X) leaves its ψ-sites occupied and its χ-sites
     empty. A left-pointing chord does the opposite.
   - Each chord therefore contributes e^{∓(p_ψμ_ψ − p_χμ_χ)/2}: only μ_eff = p_ψμ_ψ − p_χμ_χ couples to the net
     orientation m = n_→ − n_←.
   - Uncovered sites give cosh(μ_ψ/2)^{N−kp_ψ} cosh(μ_χ/2)^{N−kp_χ}. Doubly covered sites are subleading (BBNR's case 3).
4. *Fixed charge.* Fourier transforming in (iμ_ψ, iμ_χ), the two Gaussian integrals give
   e^{−2(A²+B²)/N} · e^{λ s m} · e^{−λm²/4}.
   - The first factor is the sector-dimension fraction.
   - The last factor cancels BBNR's q^{−m²/4} (4.14), because the λ there is the same combined λ.
   - What remains is BBNR (4.24)–(4.27) with charge s. ∎

*Matter.*
- Consider an operator O = Σ η_{I'K'} Ψ_{I'}†X_{K'} with independent Gaussian η, |I'| = p'_ψ, |K'| = p'_χ.
- An H-chord crossing it carries overlaps Poisson(p_ψp'_ψ/N) and Poisson(p_χp'_χ/N). So the crossing weight is
  q̃ = e^{−2(p_ψp'_ψ + p_χp'_χ)/N} = q^δ, with δ = (p_ψp'_ψ + p_χp'_χ)/(p_ψ² + p_χ²).
- O changes (A, B) by (p'_ψ, −p'_χ), i.e. s by the same δ. It is therefore a *chiral* single-chord operator in the sense
  of D7 (charge shift = crossing exponent).
- Uplift operators (D8 §1): E = ψ†χ has δ = 3/5, K = ψ†ψ† has δ = 4/5, Z = ψ† has δ = 2/5. Their daggers lower s by the same
  amounts (mirror s → −s, D7's g family).

*Exact chirality.*
- Q₀, E, K, Z are polynomials in the mutually anticommuting set {ψ_i†, χ_k}. Hence [Q₀,E] = [Q₀,K] = [E,K] = 0 and
  {Q₀,Z} = 0 identically (checked to 1e−15).
- For BPS h, O h is therefore Q₀-closed: it lies in BPS ⊕ exact, the B family. O†h lies in BPS ⊕ coexact, the A family.
  This is the support the chord measure assigns.
- *Clarification of D8 §1:* the identity [E,K] = {Q₀,Z} derived there from nilpotency holds as 0 = 0. No D8 result needed
  it to be non-trivial; {Q₀, 𝒵} = E P_H K − K P_H E (D8 §9) is unaffected.

**1.2 Symmetrized couplings: twisted chords.**

*Covariance.* For C = (G + G^{T_jk})/2 (E|G|² = 2) the antisymmetrized couplings A_ijk = C_ijk − C_jik have
  ⟨A_ijk Ā_lmn⟩ = (δ_ilδ_jm − δ_imδ_jl)δ_kn + ½[δ_ilδ_jnδ_km − δ_imδ_jnδ_kl − δ_jlδ_inδ_km + δ_jmδ_inδ_kl].
- The first group is the generic (untwisted) chord.
- The second pairs a Q-term with a Q†-term in which a ψ-index and the χ-index are exchanged: a **twisted chord**.

*What a twisted chord does.* Written out, ψ_i†ψ_j†χ_k ··· χ_j†ψ_kψ_i is a number factor on site i times ψ_j†χ_j† times
χ_kψ_k. So it moves a ψχ pair from site k to site j.

*Counting.*
- In a trace, or a diagonal matrix element in a Fock state, every mode must balance. The pair hops of the twisted chords
  in a diagram must therefore compose to the identity.
- For t twisted chords this imposes at least t extra index identifications (j = k for t = 1). Each costs 1/N relative to
  an untwisted diagram.
- **Twisted chords are therefore O(1/N) in every trace-level ("black-hole") quantity:** the density of states, the BPS
  fraction, and BPS-averaged correlators of operators with independent couplings. At leading order the symmetrized and
  generic models have the same chord diagrams.
- This is the chord-level version of CCSY §3.3, now covering all moments and correlators rather than the saddle alone.
- [num] It is consistent with: identical BPS dimensions at N' ≤ 8 (D8 §3); identical L1/L2 laws (D8 §6); the same
  dictionary accuracy in both models; and identical black-hole E–K emission coherence in both models (§5, task `twisted`).

*On graviton states.* On V-descendants V^n|Ω⟩ = Σ_{|S|=n} Π_{i∈S}ψ_i†χ_i†|Ω⟩, a pair hop maps tower configurations to
tower configurations with equal amplitude. Twisted chords then act coherently (Σ over N² hops, no identification): O(1).

*The uplift operators.*
- In the symmetrized model K_ij = E_ji − E_ij (D8 §1), so E[Ē_{i'k'}K_ij] = v(δ_{i'j}δ_{k'i} − δ_{i'i}δ_{k'j}).
- The E–K cross term of a two-point function is a twisted matter chord. Its net effect on the state is the pair creation
  ψ_i†χ_i†, a V insertion.
- It can be O(1) only between components related coherently by V, e.g. h₀₀ ∝ V h₁₁.

### 2. The black-hole window [derived + num]
- **The window.**
  - D(s) (BLY 3.17) is the zero-length-wormhole probability for |s| < 1/2 and vanishes at s = ±1/2. Outside the window
    there are no ground states: the exact count (BBNR 5.25) is nonzero only for |s| < 1/2.
  - Q shifts s by 1, so every Q-complex has at most one window sector, or two at exact edges s = ±1/2, where BBNR's count is
    ambiguous. This is CCSY's R-charge concentration in double-scaled form.
  - In charge space the window is the band |2(N_ψ − N/2) − (N_χ − N/2)| < 5/2: slope dN_χ/dN_ψ = 2 through the center.
- **Relation to CCSY's line [num, evaluated here from CCSY 3.14].**
  - CCSY's large-N line at p = (2,1) passes through the center with slope dn_χ/dj = 0.392, against 2/5 for the band.
  - The band is therefore the double-scaled linearization of the true line. The 2 % slope difference is a finite-p effect.
  - At N ≤ 10 it shifts the band by ≲ 0.2 sector, except near the boundary of charge space.
- **The tower relative to the window.**
  - The tower sector (n,n) has s = (n − N/2)/5: inside the window iff |n − N/2| < 5/2.
  - Following N at fixed n, the tower crosses the window, from s > 1/2 to s < −1/2, during about 10 uplifts. It then stays
    BPS outside, where the chord model has no BPS states.
- **Gapless edges [derived].**
  - By BBNR (5.3), with D7's sign, the exact (B) continuum at s starts at 2q^{−1/2}[cosh λ(s − ½) − 1] and the coexact (A)
    continuum at 2q^{−1/2}[cosh λ(s + ½) − 1], in units of κ.
  - So the B continuum reaches zero energy at s = +½, the A continuum at s = −½.
  - Near the edge the B-family density in φ tends to a constant while Λ ∝ φ². The edge is therefore a hard edge with
    ρ(E) ∝ E^{−1/2}.
  - The two edges are the two halves of one multiplet. Q maps the coexact family at s = −½ onto the exact family at s = +½
    of the same complex, and the two share their nonzero spectrum. This is the gapless Q = ½ supermultiplet of the N=2
    super-Schwarzian (HITZ App. A, eq. A.14, Fig. 5): it exists when the charge lattice contains ½. For p = (2,1) that is
    odd N (s ∈ (1/5)ℤ + 1/10). Example: at N = 7, (2,3) → (4,2) carries the near-zero level 0.003175 in both sectors.
- [num] *Results:* §5.

### 3. The uplift as emission into the window [derived + num]
**Bookkeeping.** Let a new BPS state sit in sector (A,B) of the N'-model, at chord charge s' = s_{N'}(A,B).
- *Slots.* For p = (2,1) the old slots sit at s' + 1/10 (00), s' − 1/10 (11), s' − 3/10 (10) and s' + 3/10 (01).
- *Dressing terms* (D8 §7):
  - the currents J₁₀ (E h₀₀, K h₁₁, Z h₀₁) land in old (A+1, B−1) at s' + 7/10;
  - J₀₁ (K†h₀₀, E†h₁₁, Z†h₁₀) lands in (A−2, B) at s' − 7/10;
  - z₁₁ (E h₀₁) at s' + 9/10, x₀₀ (E†h₁₀) at s' − 9/10;
  - z₀₀ (K h₀₁) at s' + 11/10, x₁₁ (K†h₁₀) at s' − 11/10.

**Emission spectrum [derived from D7].** Each term's normalized spectral measure is the Al-Salam–Chihara measure with
a = q^{s_in − ½} and b = q^{½ + δ − s_out}.
- *Point mass.* For s_in < ½ it has a point mass at zero energy: the BPS→BPS (obstruction) weight
  ‖P_H O h‖²/‖O h‖² = (q^{1−2s_in};q²)/(q^{1−2s_out};q²) (BLY 4.14). This is the obstruction of D8's Theorem A, now with a
  chord value.
- *Continuum.* The continuous part integrated against 1/Λ is the dressing building block, in units c₀/κ, where:
  - c₀ is the exact τ = 0 value;
  - κ comes from the exact first moment of the target sector, as in D7.

**Resonances [derived].**
- When a term lands on its gapless edge (raising at s_in = +½, lowering at −½), f(φ) → const while Λ ∝ φ². The time
  integral ∫ f/Λ dφ ∼ ∫ dφ/φ² then diverges.
- The resonant new charges are s' = −1/5 (J₁₀), −2/5 (z₁₁), −3/5 (z₀₀), and their mirrors.
- For p = (2,1), s' ∈ (1/5)ℤ at even N' and (1/5)ℤ + 1/10 at odd N'. So the currents are exactly resonant only at even N'.
  At odd N' every bulk sector is 1/10 from an edge (near-resonant).
- Prediction: black-hole dressing grows as the emission approaches a gapless edge. At finite N it is controlled by the
  near-zero modes of the edge sector.

[num] *Results:* §5.

### 4. The tower: no emission [derived + num]
- **Mechanism (D8 §4).** The Ward identity [E,V] = K and EΩ = 0 give E V^nΩ = nK V^{n−1}Ω, so J₁₀ = J₀₁ = 0. The tower emits
  nothing, at any coupling and any N.
- **Chord reading.**
  - The cancellation is between the E-emission of the 00 slot and the K-emission of the 11 slot. Their interference term is
    the twisted E–K matter chord, i.e. a V insertion (§1.2).
  - It cancels the direct terms only along a V-chain that ends on an E-annihilated state (here Ω).
  - Black-hole states are not V-descendants of E-annihilated states. On them twisted chords are O(1/N), so their
    emissions add incoherently.
- [num] *Results:* §5.

### 5. Results [num] (filled from `results/data/two_flavor_chords_2026-10-25/summary.json`)
**Window (task `window`; N' = 4–8 every sector, both models, seeds 1–3 (N' = 8: 1–2); N' = 9, 10 sectors of dim ≤ 3000).**
- *Support.*
  - At N' = 4–8, every sector with BPS states lies in |s| ≤ 1/2, except the diagonal (tower) sectors with
    |n − N'/2| ≥ 3.
  - At N' = 8 these are (0,0), (8,8), (1,1) and (7,7), the latter two at s = ∓0.6. Each of (1,1) and (7,7) holds 8 BPS
    states: the tower plus 7 forced by dimension counting (dim (1,1) = 64 > dim (3,0) = 56). The in-window partners of
    these complexes, (3,0) and (5,8), are then empty.
- *The tower alone (N' = 9, 10).* (1,1) has exactly one BPS state, the tower, with symmetrized couplings and none with
  generic ones (seeds 1–2; dimension counting gives 81 < 84 and 100 < 120). This is the first sector where the BPS
  dimensions of the two models differ; D8 found them equal for N' ≤ 8.
- *Boundary miss (N' = 10).* The 15 states of the complex J = 6 sit in (2,2), at s = −0.6, in both models. This is the
  curvature of the true line near the boundary of charge space (the slope at the center is 0.392 vs 2/5).
- *Fractions.* Exact/D(s) over window sectors of dim ≥ 20 averages 0.91–1.04 (range 0.79–1.76), and 0.84–0.95 in central
  sectors.
- *Continuum edges.* In bulk sectors (dim ≥ 50, off the edge) the exact lowest eigenvalue of each family, divided by the
  chord edge (§2), has geometric mean 0.65–0.97 and log-scatter 0.2–0.47, over gaps spanning three decades. The exact gaps
  are systematically lower, as in D7.
- *Gapless edges.* Edge sectors (s = ±1/2) of large dimension have the predicted family down to ~5×10⁻⁵ of the sector
  mean at N = 7, against a bulk median of 0.06–0.12.

**Edge-gap scaling (task `edgegap`; largest edge sector, both models).** Relative gap ε/m1 of the gapless family:
- N = 3 (dim 3): 0.05–0.27;
- N = 5 (dim 50, family 48): 4×10⁻⁴ – 4×10⁻³;
- N = 7 (dim 735, family 686): 2.9×10⁻⁵ – 1.0×10⁻⁴;
- N = 9 (dim 10584, family 9297): 3.3×10⁻⁵ – 4.8×10⁻⁵ (both models, seeds 1–2; symmetrized 3.6–4.8, generic 3.3–3.6).

Up to N = 7 the gap falls by about two orders of magnitude per step, faster than 1/D. At N = 9 it **stalls**, while the
levels above it pile up: the six lowest lie within a factor 1.4–1.55 of each other (N = 7: a factor 5). So the chord model's
exact gaplessness is approached but not reached at p = (2,1).

[conj] A possible reading: the true gapless edge (CCSY's curved line in IR-charge terms) sits at s = ½ + Δs. The chord
edge then leaves a residual gap ≈ λ²Δs²/(1 + q) in units of the mean. Δs ≈ 0.006 reproduces N = 9 and is consistent
with N = 7; at N ≤ 5 the discreteness of the levels dominates. A fixed Δs makes the residual gap fall only like
λ² ∝ N⁻² (×0.6 from N = 7 to 9), as against the exponentially fast fall of the level spacing up to N = 7. Δs is fitted,
and two points cannot fix a scaling. Whether the stall is a finite-p offset or a genuine finite-N gap needs larger N
(sparse methods).

**Dictionary (task `blocks`; N = 3–7, both models, seeds 1–3 (N = 7: 1–2), η-averaged with random phases).** Chord/exact,
geometric mean and log-scatter:
- obstruction weights (point masses): 1.05–1.21, 0.06–0.29;
- time integrals with targets ≥ 0.15 from an edge: 0.80–1.08, 0.17–0.34 (a mild drift to 0.80–0.88 at N = 7, with the
  exact gaps lying below the chord edges);
- targets 0.1 from an edge: 0.89–0.96, 0.27–0.71.

Both models agree equally well, as §1.2 requires. Targets exactly on a large edge (chord: ∞) give exact T/c₀ medians of
0.46–0.53 (N = 5) and 1.2–2.8 (N = 7), against 0.011–0.062 elsewhere. So they are enhanced 25–250×, and the enhancement
grows with N.

**Resonances (task `resonance`; exact bare fractions μ, κ = 1/μ − 1; N' = 4–7 seeds 1–3, N' = 8 seeds 1–2).**
- *N' = 6, bulk sectors.*
  - J-resonant (s' = ±1/5): the median over sectors of the mean dressing is 4.8 (symmetrized) and 9.8 (generic), against
    0.42 and 0.51 for non-resonant sectors.
  - Median per-sector maxima: 270 / 778 against 6.7 / 7.9.
  - z₁₁-resonant (±2/5): 1.2 / 1.6.
- *N' = 8 (seeds 1–2), bulk sectors, median over sectors (sym / generic):*
  - J-resonant with a large edge target: mean κ 125 / 48, largest 7.8×10⁴ / 1.5×10⁴;
  - J-resonant with a small boundary edge target (same s', no near-zero modes): 1.3 / 1.2 and 35 / 31, i.e. no
    enhancement (a control);
  - z-resonant: 2.4 / 2.2 and 34 / 28;
  - non-resonant: 1.0 / 0.98 and 30 / 29.
  So the J-resonance raises the mean ×50–125 and the maximum ×10³–10⁴, and it grows by one to two orders of magnitude
  from N' = 6.
- *Odd N' (near-resonant):* intermediate (N' = 7: mean 1.8, maxima ~44).
- Medians over states move little: the resonance lives in the upper tail.
- *D8 connection.* This explains D8 §9.4's smallest σ_min(𝕄) = 0.013 at N' = 8 in (3,3) and (5,5), with "old gap 0.0032".
  Those are the J-resonant sectors, and 0.003175 is the near-zero mode of the edge sector (4,2) at N = 7.

**Tower (task `tower`; symmetrized, N = 3–7).**
- Every current cancels: |J|/|one term| ≤ 2.5×10⁻¹⁶ over 96+ records.
- G-weighted, the individual emissions are of the size of the black-hole averages of the same operators (along the
  diagonal within a factor ~3).
- At the gapless targets they share the enhancement. At N = 7, seed 1, the E- and K-terms into (4,2) are 16.9 and 28.2 per
  unit norm, against black-hole means of 61.5 and 21.5. Off-edge they are 0.02–0.9.

**Twisted chord on black holes (task `twisted`; N = 4–7, seeds 1–3).**
- The quantity is the mean squared G-cosine between the E-emission of B(n,n) and the K-emission of B(n−1,n−1), with the
  tower projected out.
- It is the same in the two models at every (N, n). For example 0.0064 vs 0.0065 (N = 7, n = 4) and 0.0219 vs 0.0210
  (N = 6, n = 4).
- In the bulk it falls with N: 0.044 (N = 5, n = 3), 0.022–0.031 (N = 6), 0.006–0.010 (N = 7, n = 4, 5).
- The tower's own cosine² is 1 to 1e−15.

**Pairing (task `pairing`).**
- The sector means of c = ⟨n_ψe n_χf⟩ − ⟨n_ψe⟩⟨n_χf⟩ over fortuitous states lie in [−0.023, 0.021] for N' = 4–7
  ([−0.002, 0.008] at N' = 8). The symmetrized model's are slightly larger (median 0.008 at N' = 7, 0.002 at N' = 8,
  against 0.003 and 0.0007 generic): an O(1/N) twisted-chord effect, two orders of magnitude below the tower.
- Per-state 5–95 % ranges are within ±0.1, shrinking with N, in both models.
- ⟨n_ψe⟩ = A/N' to ≤ 0.035 in bulk sectors.
- The tower has c = n(N' − n)/N'² = 0.11–0.25 exactly.

### 6. Monotone vs fortuitous in the chord picture [summary]
**Does the distinction emerge in a purely chord picture? Largely yes, with one gap.**
- *From the chord rules alone*, tested in both models:
  - where black holes are: the window and its gapless edges;
  - how they lift: emission of chiral matter chords, with ASC weights and edge resonances;
  - that symmetrization does nothing to them: twisted chords are O(1/N);
  - that they absorb the new pair flavor by flavor.
- *The tower is the exception to each statement.* It sits outside the window; it emits nothing, although each emission
  alone is black-hole-sized; twisted chords are O(1) for it and supply exactly the cancelling term; and it absorbs the
  pair bound.
- *Not yet chord-intrinsic:* the reason for the exception. Zero emission is the microscopic V Ward identity. The chord
  picture identifies the cancelling term as a twisted chord (a V insertion) and shows why black holes cannot use it, but
  does not derive the identity from chord dynamics. That needs a chord Hilbert space for graviton states.

| | fortuitous (black hole) | monotone (tower, graviton gas) |
|---|---|---|
| chord object | zero-length wormhole, \|s\| < 1/2 | V-descendant of Ω; outside the window: alone, gapped |
| where | band around CCSY's line (slope 2) | diagonal N_ψ = N_χ, crossing the band |
| uplift | emits chiral matter chords E, K, Z into the window continuum | emits nothing: J₁₀ = J₀₁ = 0 |
| obstruction weight | ASC point mass (BLY 4.14) | cancelled by the 11 slot |
| dressing | time-integrated ASC two-point function; diverges at gapless edges (s' = ±1/5 resonances) | 0 at every N, including resonant sectors |
| twisted chords (symmetrization) | O(1/N): same as the generic model | O(1): the cancelling V insertion |
| absorbing the new pair | independent, c ≈ 0 | bound pair, c = n(N'−n)/N'² |

### Status
- **Derived (double-scaled):**
  - the two-flavor dictionary (λ, s, δ);
  - exact chirality, and the triviality of [E,K] = {Q₀,Z};
  - twisted-chord counting;
  - the window and its gapless edges;
  - the resonance positions.
- **Numerical, at p = (2,1), N ≤ 10:**
  - window support, and isolation of the tower at N' = 9, 10;
  - dictionary accuracy;
  - edge-gap collapse and its N-dependence;
  - dressing resonances in both models;
  - the tower's individual emissions are BH-like and cancel;
  - no twisted coherence for black holes;
  - pairing.
- **Interpretation:**
  - zero-length wormholes = black holes;
  - window edges = the super-Schwarzian Q = ±1/2 supermultiplet (HITZ App. A);
  - the tower = a graviton gas, dressed only through a symmetry.
- **Open:**
  - chord evaluation of 𝒵, α, γ (the exact uplift from chords);
  - a chord Hilbert space for graviton states (twisted chords resummed);
  - per-sector prediction of the mean exact dressing, which needs the slot weights and the interference;
  - the large-N form of the resonant growth;
  - N' = 10 at resonant parity (needs sparse methods).

---

## D10. The uplift anywhere in the window: where the jump starts and where the dressing lands

*2026-10-26. Context: D6 (the uplift at half filling), D8 §2 (Theorem A, two flavors), CV §§3–4 (LES, walks),
Kanazawa–Wettig 2017 (KW, Hodge pieces per degree; `literature/notes/kanazawa_wettig_2017.md`).
- Code: `src/uplift.py` (`OldModel(F=...)`, `harmonic_basis`, `uplift_bookkeeping`, `connecting_form`, `OldModel.dyson`);
  `src/chord_dressing.py` (`chord_sectors(..., F=)`).
- Scripts: `scripts/run_uplift_window.py`, `scripts/analyze_uplift_window.py`, `scripts/plot_uplift_window_map.py`.
- Tests: `tests/test_uplift_window.py`.
- Data: `results/data/uplift_window_2026-10-26/`. Figures: `results/figures/uplift_window_2026-10-26/`.
- Note: `research/notes/uplift_window_2026-10-26.md`.*

D6 assumed half filling at one point only: "Ea is exact because degree F+2 lies outside the window". This section
removes that assumption.
- The uplift is classified by where its dressing lands.
- The obstruction (CV's connecting map) is computed exactly in every sector, and the exceptional states are identified
  as the radical of a bilinear form.
- The one genuinely new regime is the one-step-off-centre sector, where the dressing resonates with the gapless
  window-edge multiplet.

### Setup
- Old model: N modes, supercharge Q₀ = ω₀∧ (q-form, q odd). New model: N′ = N + 1 modes, sector of fermion number F,
  level ℓ′ = 2F − N′ (CV eq. 29). The window is |ℓ′| ≤ q.
- As in D6: Ψ = a + b∧e. Here a ∈ Λ^F(o) (n_e = 0) sits at old level ℓ′ + 1 (the down step of CV's walks), and
  b ∈ Λ^{F−1}(o) (n_e = 1) sits at old level ℓ′ − 1 (the up step).
- E = η∧ with η the (q−1)-form of couplings touching e. Since q − 1 is even, E commutes with Q₀ and E† with Q₀†.
  E raises the level by 2(q−1); E† lowers it by 2(q−1).
- Hodge pieces of the old model: H (harmonic = BPS), Ex (exact = KW's V⁻), Co (coexact = KW's V⁺).
- Slot spaces: H_X ⊕ Ex_X ⊕ Co_X = Λ^F(o) and H_Y ⊕ Ex_Y ⊕ Co_Y = Λ^{F−1}(o).

### 1. Level bookkeeping and the three regimes [derived]
- **h family** (sources h ∈ H_X): the first-order dressing K h = Q₀⁺E h.
  - E h has level ℓ′ + 2q − 1 (degree F + q − 1). Q₀⁺ brings it back by −2q to level ℓ′ − 1, the n_e = 1 slot.
  - The intermediate supermultiplet is coexact(ℓ′ − 1) ↔ exact(ℓ′ + 2q − 1).
- **g family** (sources g ∈ H_Y): K′g = (Q₀†)⁺E†g passes through level ℓ′ − 2q + 1 (degree F − q) and returns to
  ℓ′ + 1. The multiplet is coexact(ℓ′ − 2q + 1) ↔ exact(ℓ′ + 1).
- **Chord charges** s = ℓ/(2q): h has (s_out, s_in) = ((ℓ′+1)/2q, (ℓ′+2q−1)/2q); g is the mirror,
  ((1−ℓ′)/2q, (2q−1−ℓ′)/2q).
  - At q = 3 and half filling these are D7's (1/6, 5/6) for even N′, and (0, 2/3), (1/3, 1) for odd N′.
- **Regime of a family**, fixed by where its dressing target lies (equivalently s_in against ½):

| regime | target | s_in | h family at | g family at |
|---|---|---|---|---|
| gapped | outside the window | > ½ | ℓ′ ≥ 2 − q | ℓ′ ≤ q − 2 |
| resonant | the window edge ±q | = ½ | ℓ′ = 1 − q | ℓ′ = q − 1 |
| obstructed | inside the window | < ½ | ℓ′ = −q | ℓ′ = q |

- So the **central** sectors |ℓ′| ≤ q − 2 have both families gapped; this is the D6 regime.
- The **off-centre** sectors |ℓ′| = q − 1 have one resonant family.
- The **edge** sectors |ℓ′| = q have one obstructed family. The other slot is empty: it lies outside the window.
- For q = 3 this means central ℓ′ ∈ {0, ±1}, resonant ℓ′ = ±2 (even N′ only), and edge ℓ′ = ±3 (odd N′ only).
- Consistency: CV's necessary conditions (34)–(35) for nonzero connecting maps, ℓ ≤ 2 − q and ℓ ≥ q − 2 on the old
  class, are exactly "target inside or at the edge of the window".
- `uplift.uplift_bookkeeping`; `tests/test_uplift_window.py::test_bookkeeping_regimes`.

### 2. The uplift on admissible sources [derived; verified]
- The BPS equations of D6 §1 do not depend on the sector.
- The only half-filling input of D6 §2 was that E a is exact. In general:
  - E a is Q₀-closed, since Q₀ a = 0 and [E, Q₀] = 0.
  - E maps exact forms to exact forms: E Q₀ x = Q₀ E x.
  - So with a = h + x (x exact), P_H E a = P_H E h. Likewise P_H E† b = P_H E† g.
- The equation Q₀ b = −εE a is solvable iff P_H E h = 0, and Q₀† a = εE† b iff P_H E† g = 0.
- **Theorem (general sector).** Define the connecting maps δ_h = P_H E P_H : H_X → H^{F+q−1}(o) and
  δ_g = P_H E† P_H : H_Y → H^{F−q}(o), and the admissible sources A_h = ker δ_h, A_g = ker δ_g. Then
  $$B(N',F)=\Psi(A_h,0)\oplus\Psi(0,A_g),\qquad a=(1+K'K)^{-1}(h+\varepsilon K'g),\quad b=g-\varepsilon K a .$$
- Proof.
  - Any BPS Ψ has harmonic parts h = P_H a and g = P_H b. They must be admissible, and the rest is fixed uniquely by
    the anti-Hermitian system of D6 §5.
  - 1 + K′K is invertible in every sector. K′ maps into exact forms, and K′K|_Ex = L†L ≥ 0.
  - The map (h, g) ↦ Ψ is injective because h = P_H a and g = P_H b.
- dim B = dim ker δ_h + dim ker δ_g. Since δ_g is the adjoint of the connecting map H^{F−q} → H^{F−1}, this is CV
  eq. (26).
- D6 §5 holds verbatim with B_o^F → A_h and B_o^{F−1} → A_g. Its proofs used only [E, Q₀] = 0 and the Hodge
  decomposition. This covers:
  - the adjoint relation P_Ex K′ P_Co = L† and the anti-Hermitian system;
  - the orthogonal splitting B = B_h ⊕ B_g and the minimum-norm uplift;
  - the closed-form decoder (𝒩, 𝒢) and monotone self-dressing.
- **Verified** (seeds 9000–9002; q = 3 at N′ = 5–14, all |ℓ′| ≤ 4; q = 5 at N′ = 6–14, all |ℓ′| ≤ 6):
  - d_h + d_g = the directly computed BPS dimension of the N′-mode model = the known count (strip walks plus KW
    exceptional states);
  - uplift residual outside the BPS space ≤ 6.4·10⁻¹³, with full rank;
  - closed-form vs direct decoder ≤ 2.6·10⁻¹¹;
  - splitting overlap ≤ 1.1·10⁻¹³;
  - monotone self-dressing (ordered h-eigenvalues at or above first order, g at or below) in all 50 single-family
    sectors, to round-off (≤ 3·10⁻¹¹).
- At N′ = 15 the code reproduces D6's exact half-filling T₂ (0.5866 and 0.6013 for seeds 9000 and 9001; D6 table:
  0.587, 0.601). At N′ = 16 it gives T₂ = 0.6075 for seed 9000 at half filling, D6's pre-registered exact value.
- One more identity: per realization, spec λ(F) = 1 − spec λ(N′ − F) **[derived; verified]**.
  - The Hodge star conjugates H_C to H_{C̄} and sends n_e to 1 − n_e; complex conjugation maps H_{C̄} back to H_C.
  - Verified on 39 (q = 3) and 40 (q = 5) sector pairs, to ≤ 2.7·10⁻¹².
  - D6's per-realization symmetry at even N′ is the self-mirror case ℓ′ = 0.

### 3. The connecting maps: exact ranks, the exchange, and exceptional states as a radical [derived; verified]

**Pairing.**
- For harmonic α at level ℓ and β at level −ℓ − 2(q−1) of the old model, define B_η(β, α) = ∫β∧η∧α (top-form
  coefficient; `uplift.connecting_form`).
- β ↦ conj(⋆β) is an antilinear bijection between harmonic spaces at levels ℓ and −ℓ (by the previous paragraph). So
  the matrix of δ on H(N, ℓ) is the matrix of B_η.
- The two ranks agree: rank δ|H(ℓ) = rank δ|H(−ℓ−2(q−1)).
- Moving β through α∧η costs (−1)^{k}, where k is the degree of α (η has even degree). On the self-paired level
  ℓ = 1 − q (k = (N − q + 1)/2):
  $$B_\eta(\alpha,\beta)=(-1)^k B_\eta(\beta,\alpha).$$

**Exact ranks for q = 3** (e(N) = exceptional states of the N-mode model per edge):

| ℓ′ | δ_h : source → target | rank δ_h | A_h | A_g |
|---|---|---|---|---|
| 0, ±1 | target outside window | 0 | H_X | H_Y |
| −2 | H(N,−1) → H(N,+3) | e(N): surjective | codim e(N) | H(N,−3): the e(N) exceptional states |
| −3 | H(N,−2) → H(N,+2) | 3^{N/2−1} − e(N′) | radical of B_η | 0 |
| −4 | H(N,−3) → H(N,+1) | e(N): injective | 0 | 0 |

- The positive ℓ′ are the mirror, with δ_g.
- Derivation:
  - The ℓ′ = −4 row is forced by the window premise at N′: no BPS states at |ℓ′| > q (KW 6.15; CV eq. 12).
  - The pairing B_η makes H(N,−3) and H(N,−1) dual, so the ℓ′ = −2 rank equals the ℓ′ = −4 rank. This gives the
    surjectivity.
  - The ℓ′ = −3 row is the self-paired form.
- **The exchange at ℓ′ = ∓2 (N ≡ 1 mod 4) [derived].**
  - e ordinary classes of H(N, −1) are obstructed: they would have to dress into the exceptional states at +3.
  - The e exceptional states at −3, which cannot step down (row ℓ′ = −4), step up instead and replace them.
  - The total 3^{(N−1)/2} is unchanged, which is CV's walk count.
  - So CV's strip walks count the dimension correctly, but an exceptional state's history is a boundary excursion.
    It is created at the edge by the radical (below), and it leaves by an up step that displaces one obstructed
    ordinary class.
- [verified] All connecting-map ranks equal these predictions in all 142 q = 3 records (N′ ≤ 16) and all 114 q = 5
  records (with the q = 5 analogues: e(N) at |ℓ′| = 4, 6; radical at |ℓ′| = 5).

**Exceptional states = the radical of an (anti)symmetric form [derived; verified].**
- The edge sector ℓ′ = −q has empty b-slot (level −q−1). So B(N′, −q) = ker δ_h on H(N, 1−q), which is the radical
  of B_η restricted to the self-paired level.
- With k = (N′ − q)/2, B_η is antisymmetric iff k is odd. An antisymmetric form on an odd-dimensional space has a
  nonzero radical. Hence:
  $$k=\tfrac{N'-q}{2}\ \text{odd}\ \ \text{and}\ \ \dim H(N'-1,1-q)\ \text{odd}\ \Longrightarrow\ \text{exceptional states at }N'.$$
- For q = 3, dim H(N, −2) = 3^{N/2−1} is always odd, so the condition is N′ ≡ 1 (mod 4). This is the uplift version
  of KW's Kramers argument (§6.2) and of CV eq. (15).
- Parity equivalence with CV's version: CV's M on Λ^a(N′) with a = k gives dim harm ≡ D̃_a(N′) (mod 2). Ours gives
  S_{N′−1}(1−q) = D̃_a(N′−1) − D̃_{a−1}(N′−1), while D̃_a(N′) = D̃_a(N′−1) + D̃_{a−1}(N′−1). These have the same
  parity.
- [verified, `forms_p*.jsonl`] |B^T − (−1)^k B|/|B| ≤ 4.3·10⁻¹⁵ for q = 3 (N = 4–14) and q = 5 (N = 6–14).
  - q = 3: radicals 1, 3, 1 at N = 4, 8, 12, and 0 at N = 6, 10, 14.
  - q = 5: radical 1 at N = 6, 10. At N = 14 the form is antisymmetric (k = 5), but dim = 1000 is even, so the radical
    is 0, as KW report for N_c = 15.
  - The singular values of the antisymmetric forms come in equal pairs, as they must.
- **The parity criterion beyond KW's range [derived].** Forced exceptional states for N′ ≤ 45:
  - q = 3: every N′ ≡ 1 (mod 4).
  - q = 5: N′ = 7, 11, 19, 23, 31, 35, 43.
  - q = 7: 9, 21, 25, 29, 37.
  - q = 9: 11, 15, 23, 39, 43.
- KW checked q̂ = 5 only up to N_c = 17 ("saturated only for N_c = 7 and 11"). **N′ = 19 at q = 5 is rigorous for
  generic couplings:**
  - CV's M on Λ⁷(19) is antisymmetric on the even dimension 50388;
  - exact(7) = rank(Q|Λ²) = 171, verified for two realizations (σ_min/σ_max = 0.89), so generically;
  - hence dim H⁷ ≡ 0 − 171 ≡ 1 (mod 2).
- **Not explained:** why the radical is larger than parity requires at N′ = 9, q = 3 (dimension 3, KW footnote 12).
  The radical subspace changes with η, but its dimension is 3 for every η tried.

### 4. Dressing in the three regimes: the window edge resonates [derived (chords, identities); num]

**Two-point function.**
- D6 §3 holds on admissible sources with H₀⁻¹ → H₀⁺. The first-order dressing cost is the time-integrated
  BPS → continuum two-point function of E at charges (s_out, s_in).
- For a non-admissible source the τ-integral grows linearly. Its coefficient ‖δh‖² is the late-time plateau, i.e. the
  BPS → BPS weight (D7's point mass, = BLY 4.14). **The obstruction is the plateau of the two-point function.**

**Late-time behaviour in the chord limit** [derived from D7 §1]. The ASC weight is
f = (ab; q²)(q², e^{±2iφ}; q²)/[2π(ae^{±iφ}, be^{±iφ}; q²)] with a = q^{s_in−½}, and
Λ(φ) = 2q^{−½}[cosh λ(s_in−½) − cos φ].
- s_in > ½ (a < 1): f ∝ φ² as φ → 0. So G(τ) ~ τ^{−3/2} e^{−κΛ_min τ}, with Λ_min = 2q^{−½}[cosh λ(s_in−½) − 1].
- s_in = ½ (a = 1): the zero of (e^{±2iφ}) cancels against (e^{±iφ}), so f(0) is finite. Λ ≈ q^{−½}φ², so
  **G(τ) ~ τ^{−1/2}** and the time integral diverges.
  - D7's closed form ∝ ₂φ₁(q², q²; q^{2+2δ}; q², q^{2s_in−1}) has terms that tend to a constant, so it has a simple
    pole ∝ 1/(λ(2s_in − 1)) as s_in → ½⁺ at fixed λ.
  - The Schwarzian limit, (2s_in−1)^{δ−1} (D7 §2), diverges too.
- s_in < ½: a zero-energy point mass, i.e. the obstruction.

**The resonant jump law** [derived under D7's chord rules].
- Near φ = 0, g(φ) ∝ A(φ)/Λ(φ) ≈ G₀/φ², and the B-family density has a finite limit ρ₀.
- Each inner mode with a large jump contributes one large eigenvalue of 𝒟 (rank one). Hence
  $$P(\sigma^2>K)\simeq r\,w(g>K)\simeq r\rho_0\sqrt{G_0}\,K^{-1/2}.$$
- The mean is therefore infinite in the chord limit. At finite N it is cut off by the smallest inner energy.
- Parameter-free amplitude (`chord_dressing.chord_jump_law(F=(N′−2)/2)`, `chord_tail.json`): P·√K = 0.50, 0.51, 0.52
  at N′ = 12, 14, 16.

**Finite N [num]** (`summary.json`, seeds 9000–9002, `fig_resonance`):
- **Which multiplet.** The resonant family dresses through coexact(−q) ↔ exact(+q). This is KW's window-straddling
  multiplet, GOE/GSE (KW 6.13, 6.24), and Turiaci–Witten's CT-invariant multiplet (−q̂/2, q̂/2), which they find β = 1
  at N = 11 (TW Fig. 4).
  - [heuristic] It is the only multiplet in which Q₀ maps between ambient spaces of equal dimension,
    C(N, (N−q)/2) = C(N, (N+q)/2). That suggests a "square" hard edge at zero.
- **Its gap** (seed means). E_min/Ē = 2.5e−1, 2.1e−2, 2.9e−3, 5.9e−6, 1.6e−6, 8.0e−8 at N = 5, 7, 9, 11, 13, 15
  (N = 15: seed 9000; seed 9001: 1.5e−10). Gapped families sit at 0.02–0.6.
  - n²E_min/Ē, with n the multiplet dimension (4, 21, 80, 319, 1208, 4551), stays at 0.003–24 per seed (q = 3,
    N = 5–15; 1.65 and 0.003 at N = 15), and at 0.15–13 for q = 5.
  - The smallest values (0.003 and 0.06) occur in the GOE multiplets (N ≡ 3 mod 4). This fits the softer β = 1 hard
    edge, P(nσ_min < x) ∝ x. The Kramers (GSE) multiplets never go below 1.1. With 2–3 seeds per N this is
    suggestive only. So E_min ~ Ē/n² up to an O(1–10) factor, exponentially small in N.
  - [heuristic] This is the hard edge of a square random matrix: Edelman's square Gaussian law gives n²x_min = O(1).
    The prefactor fluctuates by more than a decade across seeds, with no clean GOE/GSE separation.
  - The h family at ℓ′ = −1 (s_in = 2/3) closes much more slowly: 0.21 → 0.017 from N = 6 to 14. That is a power
    of N, in line with a chord gap Λ_min ∝ λ²(s_in − ½)² at small λ = 2q²/N; the matching was not checked.
- **Dominance of the bottom modes.** The lowest 1 % of inner modes carry 90 %, 96 % and 99 % of Tr(first-order
  dressing) at N = 11, 13, 15. In gapped families they carry ≤ 14 %.
- **Dressing costs (first order; seed means).** Mean 1.3, 4.1, 7.0, 1244, 851, 1579 / 6.6·10⁵ (per seed; the latter from one mode with E/Ē = 1.5·10⁻¹⁰, and 613 / 688 without the top state) at N′ = 6, 8, 10, 12, 14, 16;
  median 0.14–0.74; maximum 1.45·10⁹ (N′ = 16, seed 9001). Central sectors: mean ≤ 0.42, maximum ≤ 5.2. The mean is dominated by single modes
  and fluctuates strongly between seeds, as a K^{−1/2} tail requires.
- **The tail.**
  - P(σ² > K) ∝ K^{−α} over four decades. Hill α (K > 100, h family at ℓ′ = −2, all seeds) is 0.43 ± 0.06
    (N′ = 12, n = 48), 0.48 ± 0.04 (N′ = 14, n = 144) and 0.48 ± 0.03 (n = 259) (N′ = 16). The chord value is α = ½.
  - Amplitude P·√K at K = 100: 0.66, 0.66, 0.59 (data) vs 0.50, 0.51, 0.52 (chord) at N′ = 12, 14, 16. The ratio
    data/chord falls with N′ (1.32, 1.29, 1.13), i.e. toward the chord law.
  - q = 5 shows the same at |ℓ′| = 4 only: gap 8.1e−5 (N = 11) and 7.5e−7 (N = 13); mean dressing 98 and 1070.

**What the exact uplift does with the resonance** [derived identity + num].
- In the singular basis Q₀ = Σ s_j |u_j⟩⟨v_j| of the inner multiplet, with e = U†E V_h:
  $$K_H^\dagger(1+R)^{-1}K_H=e^\dagger(S^2+W)^{-1}e,\qquad W=U^\dagger E P_{\rm Ex}E^\dagger U = S R S .$$
  - This is the first-order dressing with the propagator of H₀ + E P_Ex E† (a Dyson resummation; `OldModel.dyson`,
    verified ≤ 4·10⁻¹⁵).
  - But [E, Q₀] = 0 forces W_jj = s_j² R_jj with R_jj = O(0.1) (measured 0.06–0.14 on the lowest modes). The
    self-energy vanishes on the near-zero modes, so **the resonance survives in the exact (minimum-norm) dressing
    cost**: mean 1190 vs 1244 first order at N′ = 12, 787 vs 851 at N′ = 14, and 1454 vs 1579 (seed 9000) and 6.1·10⁵ vs 6.6·10⁵ (seed 9001) at N′ = 16.
- **The decoder does not resonate.**
  - For a cloud-dominated state, the exact h-block gives λ → 1 − ⟨k̂,(1+R)^{−2}k̂⟩/⟨k̂,(1+R)^{−1}k̂⟩. This is set by
    the bounded operator R, not by σ.
  - Measured: first order puts 2.0–2.5 % of resonant-sector states within 10⁻³ of the wall (lowest λ₁ ≈ 10⁻⁷–10⁻⁵).
    The exact decoder has none, and none within 10⁻² at N′ = 16.
  - The exact floor rises with N′: the lowest λ is 4.5–6.5·10⁻³ (N′ = 12), 8–10·10⁻³ (N′ = 14) and (1.4–1.5)·10⁻² (first order: 5·10⁻⁷, 7·10⁻¹⁰) (N′ = 16).
  - D6's per-mode mean field (R → ρ_j) fails here (T₂ low by 0.04–0.05), while it is good in central sectors. The
    spread of R along the resonant direction matters.
- **[interpretation]**
  - In resonant sectors a few BPS states of N′ are dominated by the near-zero modes of the old edge multiplet: the
    "almost-BPS" states of N at the window edge s = ±½.
  - Their harmonic seed carries a vanishing fraction of the norm, ‖h‖²/‖Ψ‖² ~ 1/κ.
  - This is CV's metric fragility in its strongest form, and it is localized at |ℓ′| = q − 1.
  - In chord language: the dressing two-point function of these sectors ends on the gapless Q = ½ supermultiplet of
    the super-Schwarzian (HITZ App. A).
  - Two-flavor analogue: the resonances at s′ = ±1/5 (D9).

### 5. The decoder across the window [identities; num]
Seed-averaged values (3 seeds; N′ = 16: 1–2 seeds), from `summary.json`:

| N′ | ℓ′ | d (h+g) | T₂ first → self → exact | H first / exact | T₄ − (T₂ − H_first) | walls(10⁻³) first / exact |
|---|---|---|---|---|---|---|
| 12 | 0 | 486 (243+243) | 0.631 → 0.670 → 0.685 | 0.085 / 0.090 | −0.0045 | 0.329 / 0.329 (atoms) |
| 12 | −2 | 243 (243+0) | 0.515 → 0.497 → 0.497 | 0.125 / 0.132 | −0.0069 | 0.025 / 0.000 |
| 13 | −1 | 729 (486+243) | 0.563 → 0.604 → 0.616 | 0.111 / 0.109 | +0.0026 | 0.112 / 0.123 |
| 13 | −3 | 1 (1+0) | λ = 0.65, 0.58, 0.65 (one exceptional state) | – | – | – |
| 14 | 0 | 1458 (729+729) | 0.566 → 0.619 → 0.636 | 0.105 / 0.108 | −0.0034 | 0.116 / 0.128 |
| 14 | −2 | 729 (728+1) | 0.477 → 0.455 → 0.455 | 0.133 / 0.142 | −0.0092 | 0.023 / 0.000 |
| 16 | 0 | 4374 (2187+2187) | 0.539 → 0.595 → 0.608 (1 seed) | 0.122 / 0.125 | −0.0030 | 0.011 / 0.021 |
| 16 | −2 | 2187 (2187+0) | 0.452 → 0.427 → 0.427 | 0.138 / 0.148 | −0.0106 | 0.020 / 0.000 |

- **T₄ = T₂ − H** is exact in every sector. The first-order prediction T₄ ≈ T₂ − H_first (D6 §6) holds to ≤ 0.009 in
  every extensive sector at N′ = 9–14, and to ≤ 0.011 at N′ = 16.
  - The resonant sectors are the worst, and there the error grows steadily: −0.0069, −0.0092, −0.0106 at
    N′ = 12, 14, 16 (central: −0.0045, −0.0034, −0.0030).
  - Their first-order spectrum is visibly wrong at the bottom, so this growth is expected.
  - D6's H_first ≈ H is therefore not uniform across the window: it degrades where the dressing resonates.
- **In resonant sectors self-dressing lowers T₂.** Cloud-dominated states move away from the λ = 0 wall, so y drops
  from ≈ 1. In central sectors self-dressing raises T₂.
- **Atom counting** (D6 §4) generalizes: forced atoms = max(0, dim A − dim(target family)).
  - At ℓ′ = −2 the target family is the edge multiplet. N′ = 8: 27 − 21 = 6 atoms, observed 6 in every seed. None from
    N′ = 10 on.
- **Exceptional states.** They are h-type at ℓ′ = −q (g-type at +q). The state at N′ = 5, 9, 13 has 58–76 % of its
  weight in n_e = 0 at the lower edge (spread 0.27–0.79 at N′ = 9).
- The exceptional b-states that uplift at ℓ′ = −2 have gapped first-order dressing (target −7): λ₁ = 0.02–0.13.
  - N′ = 6: an exact atom at λ = 0, since E† kills a 1-form.
  - N′ = 10: in the exact decoder they sit at λ = 0.003–0.02, with g-weight 0.54–0.94.
  - N′ = 14: the single exceptional g-state dissolves into the resonant h-block. No exact eigenvector has g-weight
    above 5–9 %. The h–g mixing block −2εK_H†L(1+R′)⁻²K′_H is driven by the large resonant K_H.

**Numerical caution [num].** The exponentially small edge gap makes rank decisions delicate.
- At N = 15 (seed 9001) the smallest edge mode (σ/σ_max ~ 10⁻⁵) fell below the Gram-eigenvalue threshold
  σ² > 10⁻¹⁰σ²_max of `bps_basis_fast`, and was counted as a harmonic form at level −3.
- The dimension check caught it (d = 2187 + 1).
- `harmonic_basis` now decides ranks by singular values (σ > 10⁻¹⁰σ_max, as `pinv` does); the record was rerun.
- Beyond N′ ≈ 18 the edge gap approaches the precision of double-precision dense linear algebra. Larger sizes need
  rank decisions informed by the known dimensions, or higher precision.

### 6. Figure (`scripts/plot_uplift_window_map.py`)
- `fig_jump_map_q3` / `_q5`: CV Fig. 2/3 conventions with N up to 16. The right margin gives the regime of each row
  ℓ′; red crosses mark sectors with δ ≠ 0.
- `fig_jump_panels_q3`: the four regimes, each drawn as one step N → N′:
  - central (11 → 12, ℓ′ = 0);
  - resonant (11 → 12, ℓ′ = −2);
  - resonant with exchange (9 → 10, ℓ′ = −2);
  - edge (12 → 13, ℓ′ = −3).
- Each panel shows the jump arrows from the two slots, the E / E† dressing loops to their targets coloured by regime,
  and the Q₀⁺ returns.

### Status
- **Exact:**
  - §2: the general-sector theorem, all of D6 §5 on admissible sources, and the mirror symmetry;
  - §3: the duality of connecting maps, the (anti)symmetry of B_η, the radical characterization of exceptional states,
    the parity criterion, and N′ = 19 at q = 5 (given rank(Q|Λ²) = 171, checked);
  - §4: the Dyson identity with W = SRS, and the cloud-dominated decoder limit;
  - T₄ = T₂ − H.
- **Exact given the window premise** (KW 6.15 / CV Conj. A.1, verified to N ≤ 17): the rank table of §3 and the
  exchange.
- **Derived under D7's chord rules:** the late-time asymptotics e^{−E₀τ}τ^{−3/2}, τ^{−1/2} and the plateau; the
  K^{−1/2} tail and its amplitude.
- **Numerical:**
  - everything in §§2–3 is verified to machine precision;
  - the edge-multiplet hard edge (E_min/Ē ≈ 10⁻⁶ at N = 11–13);
  - the tail exponent 0.43 → 0.48 (chord: ½) and amplitude within ~30 % of the chord value;
  - no wall pile-up in the exact decoder;
  - the T₄ prediction to ≤ 0.009 in all sectors.
- **Heuristic:** the "square ambient dimensions" reason for the hard edge.
- **Open:**
  - the N′ = 9 radical excess;
  - the large-N law of the resonant tail and of E_min (finite-size vs chord);
  - a closed form for the exact resonant decoder edge (λ_min ≈ 10⁻²);
  - whether the "cloud-dominated" BPS states have a gravitational description: BPS states of N′ assembled from the
    gapless edge continuum of N;
  - two-flavor and q ≥ 5 versions of the radical criterion with the full counts.

---

## D11. Chords in the resonant sectors, layer 1: the first-order dressing spectrum

*2026-10-26. Context: D7 (the first-order dressing spectrum from chords, central sectors) and D10 §4 (the window-edge
resonance). The resonant chord problem is split into three layers:
1. the first-order spectrum (this section);
2. the finite-N hard edge;
3. the self-dressing operator R, i.e. continuum → continuum chords.

Code: `src/chord_dressing.py` (`DressingFamily(inner='svd')`, `old_model_from_uplift`, `chord_population(F=)`,
`analyze_family(Fsec=)`). Scripts: `scripts/run_chord_resonance.py`, `scripts/analyze_chord_resonance.py`. Data:
`results/data/chord_resonance_2026-10-26/`. Figure: `results/figures/chord_resonance_2026-10-26/fig_layer1.pdf`. Test:
`tests/test_uplift_window.py::test_dressing_family_in_resonant_sector`.*

### 1. What the chord rules predict [derived from the chord rules]
- The resonant family is the h family at ℓ′ = −(q−1). For q = 3 its charges are (s_out, s_in, δ) = (−1/6, ½, 2/3).
- D7's rules apply unchanged:
  - uncrossed diagrams, plus BPS faces on the unique zero mode;
  - so the first-order operator 𝒟 = P E† H₀⁺ E P is free compound Poisson, with rate r = e/b and jump law
    g(φ) ∝ A(φ)/Λ(φ) under the B-family density.
- What changes at s_in = ½ is only that the jumps are unbounded: g ≈ G₀/φ² near φ = 0 (D10 §4). Consequences:
  - The moments of 𝒟 (t₁ = mean dressing, …) are infinite in the chord limit.
  - The spectrum is still a well-defined free compound Poisson law, whose large-σ² tail is
    P(σ² > K) ≈ r ρ₀ √G₀ K^{−1/2}.
  - Bounded functions of σ² are finite and insensitive to the cutoff. This includes the first-order decoder
    y = ((1−σ²)/(1+σ²))², W^(1) = ⟨y(1−y)⟩ and T₂^(1) = ⟨y⟩.
- The prediction is computed as in D7 §3. The chord population of S (n = e inner levels at quantiles of the
  B-family density) goes into the separable model X = S^{1/2} Z (factorized BPS faces, C = 1) at the exact (e, b).
  The top quantile (e − ½)/e plays the role of a finite-e cutoff on the 1/φ² tail.

### 2. Implementation notes [num]
- The inner spectrum must come from singular values. D7's Gram cut (σ² > 10⁻⁹σ²_max) would drop the near-zero edge
  modes at N = 15, the same failure as in D10 §5. `DressingFamily(inner='svd')` agrees with `'eigh'` to 3·10⁻¹⁵ at
  half filling.
- Sources are the admissible ones from `uplift.OldModel(F)`.
- The own-η spectrum equals the window runs' σ² to 2·10⁻¹⁵ (same couplings).
- Error budget as in D7 §4:
  - "sampled": the exact first order over the η ensemble;
  - "wishart": the exact inner S with C = 1, which tests the planar structure;
  - "chord": the D7 jump law, which tests the chord continuum.

### 3. Results [num] (seed means, 3 seeds; N′ = 16: 2 seeds, 8 etas)
| N′ | sector | W: sampled / wishart / chord | W error: planar, continuum, total | T₂: sampled / wishart / chord | T₂ error: planar, continuum, total |
|---|---|---|---|---|---|
| 12 | resonant ℓ′ = −2 | 0.1262 / 0.1332 / 0.1370 | +5.5 %, +2.8 %, +8.5 % | 0.5174 / 0.5020 / 0.4718 | −3.0 %, −6.0 %, −8.8 % |
| 14 | resonant ℓ′ = −2 | 0.1331 / 0.1387 / 0.1419 | +4.3 %, +2.3 %, +6.7 % | 0.4767 / 0.4650 / 0.4409 | −2.4 %, −5.2 %, −7.5 % |
| 16 | resonant ℓ′ = −2 | 0.1369 / 0.1418 / 0.1435 | +3.6 %, +1.2 %, +4.8 % | 0.4493 / 0.4370 / 0.4178 | −2.7 %, −4.4 %, −7.0 % |
| 12 | central ℓ′ = 0 (control) | 0.0853 / 0.0894 / 0.0902 | +4.9 %, +0.9 %, +5.8 % | 0.6340 / 0.6331 / 0.6143 | −0.2 %, −3.0 %, −3.1 % |
| 14 | central ℓ′ = 0 (control) | 0.1051 / 0.1101 / 0.1123 | +4.8 %, +2.0 %, +6.9 % | 0.5694 / 0.5675 / 0.5524 | −0.3 %, −2.7 %, −3.0 % |

- **The control reproduces D7** (W^(1) chord high by 4–9 %, planar about +5 %, chord inputs +1–3 %). This validates
  the generalized pipeline.
- **The heavy-tailed free compound Poisson structure holds.**
  - The planar error in the resonant sector (W +4–6 %, T₂ −2 to −3 %) is the same as in central sectors.
  - The planar model with the exact S reproduces the exact tail: P·√K at K = 100 is 0.677 / 0.650 / 0.602 (wishart)
    against 0.661 / 0.649 / 0.593 (sampled) at N′ = 12 / 14 / 16.
  - So uncrossed diagrams with factorized BPS faces still describe the spectrum when the jumps are unbounded.
- **The chord continuum is the weak link at the edge.**
  - W is as good as in central sectors (continuum +2–3 %).
  - T₂ is about twice as far off (continuum −6.0, −5.2, −4.4 % at N′ = 12, 14, 16, against −3 % central).
  - The tail amplitude ratio chord/exact is 0.74, 0.79, 0.88 (P·√K(100) = 0.49 / 0.51 / 0.52 against 0.66 / 0.65 / 0.59).
  - Exponents: the chord population gives Hill α = 0.50–0.51 (K > 100). The exact data give 0.44, 0.48, 0.48.
  - **Every error shrinks with N′.** Total W error +8.5 → +6.7 → +4.8 %, T₂ −8.8 → −7.5 → −7.0 %. The planar part of W
    falls +5.5 → +4.3 → +3.6 %, as D7 found at half filling. The chord description of the resonant first order
    converges toward the exact one.
  - The chord jump law has too little weight at small φ. In other words, the finite-N spectrum near the hard edge is
    denser than the double-scaled continuum. This is the same direction as D7's "the exact charged spectrum extends
    below the chord gap". It is the input for layer 2.
- **The Kronecker (separable) model fails in resonant sectors.**
  - It uses C = E_η𝒟/t₁, the measured BPS heterogeneity. At N′ = 12 / 14 it gives W = 0.092 / 0.101 and
    T₂ = 0.71 / 0.70, against sampled 0.126 / 0.133 and 0.52 / 0.48. In central sectors it is fine (≤ 5 %), as in D7.
  - [interpretation] E_η𝒟 = Σ_n M_n/E_n is dominated by a few edge modes with low-rank coupling patterns M_n. The
    BPS-side heterogeneity is then not separable from the inner side, and the uniform-face model (C = 1) is the right
    planar model here.
  - This is the resonant counterpart of D7 §5's coherent heterogeneity.
- **First order misses the bottom of the exact decoder** (D10 §4; fig_layer1 b). Exact and chord first order both put
  a pile of states at λ ≈ 0 that the exact uplift does not have. That is layer 3.

### Status
- **Derived (chord rules):** the heavy-tailed free compound Poisson law at s_in = ½, and finite, cutoff-insensitive
  first-order decoder statistics.
- **Numerical:**
  - the planar structure holds (tail amplitude to ≤ 2.5 %);
  - chord vs exact first order: W +8.5 → +4.8 %, T₂ −8.8 → −7.0 % for N′ = 12 → 16 (central: +6–7 %, −3 %);
  - the tail amplitude is 26 → 12 % low, a continuum effect that shrinks with N′;
  - the Kronecker model fails.
- **Next:**
  - layer 2: the finite-N hard edge, i.e. the 25 % tail deficit, E_min, and the KW class;
  - layer 3: R from continuum → continuum chords.
