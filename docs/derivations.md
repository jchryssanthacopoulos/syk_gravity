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
    5e−5 (generic) in (3,3) and (5,5) at N = 7→8, against 1e−4 to 4e−3 at N = 5, 6. Near-"pure dressing" BPS states
    exist (a lower near-wall), and injectivity may fail at larger N.
- **The E₂ page is not the answer [num].** The slot complex built from old BPS states is
  H(01) → H(00) ⊕ H(11) → H(10), with d₁ = ([K], [E]) and (h₀₀, h₁₁) ↦ [E]h₀₀ − [K]h₁₁. Then:
  - d₁² = 0, because [E][K] − [K][E] = [{Q₀,Z}] = 0, and dim B_{N+1} ≤ dim E₂;
  - equality fails in exactly four sectors per N, identically in both models.
  - At N = 6 these are (2,3), (3,5), (4,2), (5,4), with E₂ = 92, 78, 78, 92 against 14, 0, 0, 14 BPS states. At
    N = 7 they are (2,3), (4,2), (4,6), (6,5), with E₂ = 91, 330, 330, 91 against 0, 239, 239, 0.
  - A differential d₂ ∼ [Z − EQ₀⁺K + KQ₀⁺E] pairs E₂^{w=−1}(A,B) with E₂^{w=+1}(A+2,B−1), always with the maximal
    rank (e.g. 9, 78, 91 at N = 5, 6, 7).
  - So the pair uplift is not computed by the slot complex alone. The two-step formula is exact and accounts for d₂.
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

### Status
- **Exact:**
  - §1: the block form, the identities and the slot equations;
  - Theorem A, with its proof, splitting and minimum-norm uplift (§2);
  - the two-step dimension formula (§3);
  - the tower's protection, bareness and slot weights (§4);
  - the bare-map equations and the rooms (§§5–6);
  - the first-order formulas (§7).
- **Numerical laws** (N = 2–6 with 5 seeds, N = 7 with 1–2 seeds; generic control at all N):
  - four-slot injectivity (holds, but the smallest bare fraction reaches 1.8e−6 at N = 7);
  - d₂ ≠ 0 in exactly four sectors per N, with maximal rank;
  - the room laws L1 and L2, and the tower-cost criterion;
  - bare = generic count + tower;
  - BPS dimensions are coupling-independent;
  - first order within 4–6 % for sector means.
- **Unexplained:** the one-dimensional coincidences at N = 5 in the near-central sectors (2,1), (2,2), (3,3), (3,4) — two
  particle–hole pairs — in both models.
- **Open:**
  - a proof of four-slot injectivity;
  - a closed form for d₂;
  - the CV slot decoders in closed form (the analogue of D6 §5's pencil);
  - the chord computation of the two-point functions of E, K, Z and of their interference (Stage 2).
