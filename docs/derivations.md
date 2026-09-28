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

**1. Free compression is the zero-length-wormhole term [derived within DSSYK; verified numerically].** BLY give
$a = P_0$, the probability that the SUSY wormhole has chord number 0 (eq. 1.1/3.15). Our $r$ is
$\sum_n P_n\,q^{2\Delta n}$ over the same distribution (eqs. 4.5–4.6). Since $q^{0} = 1$,
$$
r \;=\; \underbrace{a}_{\text{Haar / Wachter}} \;+\; \sum_{n\ge1} P_n\,q^{2\Delta n},
\qquad
\eta_r \equiv \frac{r-a}{1-a} = \big\langle q^{2\Delta n}\big\rangle_{n\ge1}.
$$
Checked from the HH recursion (BLY 3.2–3.11) at $N'=14$: $\sum_{n\ge1}P_n = 1-D$ and
$D + \sum_{n\ge1}P_nq^{2n/3} = r_{\rm chord}$ to $10^{-16}$.
- The free-probability (Wachter/MANOVA) value of the second free cumulant is exactly the $\ell = 0$ contribution.
- The deviation from freeness, $r - a$, is the contribution of wormholes of nonzero length.
- The rigidity fraction $\eta_r$ is the average matter propagator $q^{2\Delta n}$ over those wormholes.

**[conj]** The whole free-compression law (all moments) is the $\ell=0$ truncation of the super-chord computation.
This is untested beyond $m_2$.

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
