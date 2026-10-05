# The decoder fourth moment: isotropic vs local deformations, free liberation, and the remainder R

*2026-10-05. Phase A1–A2 of the T₄/R program (plan of 2026-10-05; `research/pdfs/progress_2026-09-29.pdf` §7.3
item 2). Derivations: `docs/derivations.md` D5. Literature step: `literature/synthesis.md` §8 and the notes on BBNR,
BINT, Pluma–Speicher, Bożejko–Kümmerer–Speicher, G. Lin 2026. Code: `scripts/run_ou_fourarc.py`,
`scripts/analyze_ou_fourarc.py`, `scripts/plot_ou_fourarc.py`, tests `tests/test_ou_fourarc.py`.
Data: `results/data/ou_fourarc_2026-10-05/`. Figure: `results/figures/ou_fourarc_2026-10-05/`.*

Labels: **[derived]** proven here or in D5; **[num]** numerical observation; **[interp]** interpretation.

## 0. Summary

**Question.** The decoder fourth moment is
T₄ = τ(PP′PP′) = a² + 2a(T₂ − a) + S + R, with τ = Tr/d and P′ = UPU, U = (−1)^{n_i}. The remainder R is the part not
fixed by one- and two-copy size data. What is it, and what is its gravitational meaning?

**Answer in three layers.**
1. **Exact structure [derived, D5].**
   - The decoder, the isotropic OU family (all couplings rotated, P(ω) vs P(xω + √(1−x²)ω′)) and the site-local
     family (only couplings containing site i rotated) are **three evaluations of one generating function**,
     T(K) = Σ_α K^α t_α, at different coupling correlations K.
   - In chord language, chords between adjacent arcs carry K_J and chords between opposite arcs carry 1 (a ℤ₂ rule,
     not BINT's random-matter product rule).
   - With isotropic K this is the overlap of two mismatched two-sided BPS wormholes (an out-of-time-order
     correlator of the q-OU process).
2. **Isotropic deformations ≈ free liberation [num].**
   - Along the whole OU family, T₄ as a function of T₂ is close to the **free-liberation curve**, i.e. P rotated by a
     free unitary Brownian motion (Demni–Hamdi–Hmidi). The latter is analytic, with a closed-form R_lib(η) (D5 §5).
   - The deviation is small and of the size of the known two-model non-freeness.
   - Infinitesimally, the heterogeneity ratio κ of the principal angles is within 3–12 % of the free-Poisson value
     1/(1 − a).
   - Generic (isotropic) coupling changes thus rotate the BPS space approximately like a free unitary Brownian motion:
     a dynamical form of the BPS-chaos (Haar-orientation) hypothesis.
3. **The decoder is not isotropic: localization excess [num].**
   - At equal T₂ the decoder has substantially larger T₄, T₆ and T₈ than both OU and liberation, i.e. more
     heterogeneous principal angles.
   - The Gaussian site-local family reproduces this continuously. Its z = −1 endpoint *is* the decoder; resampling
     the site's couplings (z = 0, no ℤ₂ sign) already carries about two thirds of the excess.
   - The excess is visible infinitesimally (κ 22–49 % above free at N′ = 8–12).
   - Shares of T₄^dec − T₄^free at N′ = 8 → 14: liberation 68 → 78 %, isotropic −5 → +3 %, localization 37 → 19 %.
     The absolute localization term grows but saturates: 0.017, 0.032, 0.042, 0.047 at even N′ = 8–14.
   - So R = R_lib(η) + δ_iso + δ_loc, and δ_loc (the probe's "UV imprint") dominates the non-free part of R.

## 1. Setup and exact results (summary of D5)
- **Generating function.**
  - For a diagonal coupling correlation K, T_{2m}(K) = E τ((P(C)P(C̃))^m) with C̃ = KC + √(1−K²)C′.
  - Mehler gives T₄(K) = Σ_α K^α t_α.
  - Decoder: K_J = σ_J = (−1)^{[i∈J]} (exact: U P(C) U = P(σC), tested to 10⁻¹⁰). OU: K ≡ x. Site-local:
    K_J = z for J ∋ i and 1 otherwise, so the decoder is z = −1.
- **When decoder = OU.** The site-averaged decoder weight of a chaos multi-index is 1 − 2o(α)/N, where o is the odd
  cover of the cross-chord index sets.
  - This equals x_s^{|α|} exactly when the index sets are independent.
  - Hence decoder ≡ OU(x_s) in double scaling for heavy probes (s ∝ √N).
  - At p = 3 the index sets of overlapping ("friend") chords are strongly correlated (D4 item 6), so the matched-T₂
    comparison is the meaningful test.
- **Exact OU properties.**
  - R^OU(0) = 0 and R^OU(1) = (1 − a)² − S.
  - Sum rule Σ n θ_n = 2 Σ n π_n.
  - T₄ − 2T₂ + 1 = τ((PP_xP − P)²) = E sin⁴θ.
  - Small-angle expansion in the quantum geometric tensor: metric only, the Berry curvature cancels.
- **Free liberation (D5 §5).**
  - T₂ = a + (1−a)e^{−t}, so η = e^{−t}.
  - R_lib(η) = 2a(1−a)(1−2a)η + [E + 2a(1−a)² ln η]η².
  - This explains the observed sign changes of R at a = ½ and a = 0.618.

## 2. Numerical design
- **Coupling families** (`run_ou_fourarc.py`). Per sample, K = 16 BPS projectors are built on a circle:
  - OU: ω_j = cos θ_j ω + sin θ_j ω′ with θ_j = jπ/K, so all unordered pairs are OU pairs with x = cos(Lπ/K),
    L = 1..8. L = 8 is x = 0 exactly (independent).
  - Site-local: only the couplings containing site N′−1 rotate, θ_j = 2πj/K, giving z = cos(2πL/K) ∈ [−1, 1].
    L = 8 is the decoder.
- **Spectra.** For each pair, cos² of the principal angles = the squared singular values of V_i^†V_j, where V is an
  orthonormal BPS basis (Hodge projector + randomized range finder, ‖PV − V‖ < 4·10⁻¹⁶).
- **Matched probes.** On the same realizations, the decoder (all N′ single sites) and parity sets
  (3 random sets per |A| = 2..N′/2) are evaluated on 2 projectors per sample (1 for the site-family runs):
  μ = eigenvalues of V^†U_AV.
- **Analysis** (`analyze_ou_fourarc.py`).
  - Sample means. Bootstrap over samples (OU curve and probes resampled jointly).
  - Matched-T₂ interpolation by PCHIP in T₂.
  - Liberation by integrating the DHH ODE.
- **Seeds.** OU: 1_000_003·(seed0 + s) + 29, seed0 = 0. Site family: seed0 = 1000. Probe sets use seed + 7.
- **S.** S^OU(x) = E Σ_k φ_k(ω) w_k(ω_x) is x-independent to within the covariance of per-realization size data. Using
  `parity_family_2026-09-29/decoderT4_*.jsonl`, same-realization minus cross-realization Σφw is
  ≤ 8·10⁻⁷ for N′ = 8–14. So R^OU(x) = T₄(x) − c₀ − 2a(T₂(x) − a), with c₀ = T₄^OU(0).

## 3. Results [num]

*(Tables generated by `analyze_ou_fourarc.py --markdown --site`. Numbers in parentheses are bootstrap SEs in units
of 10⁻⁴.)*

**Figure** (`results/figures/ou_fourarc_2026-10-05/ou_fourarc.pdf`).
- (a) T₄ vs T₂ for the OU family, the site-local family and free liberation, with the decoder (★) and parity probes (▲),
  at N′ = 12, 14.
- (b) R vs η for all N′.
- (c) The three-way split of T₄^dec − T₄^free.
- (d) Principal-angle spectra at N′ = 14:
  - decoder;
  - OU at matched T₂ (histograms interpolated linearly in T₂ between neighbouring lags);
  - finite-D free liberation (D = 2000, unitary BM by 40 GUE-exponential steps, seed 11) at t = −ln η.

**3.1 Decoder (s = 1) vs OU vs free liberation at equal T₂.**
- T₄^free = 2a² − a³ (Wachter).
- Split columns:
  - liberation = T₄^lib(T₂^dec) − T₄^free;
  - isotropic = T₄^OU(x_*) − T₄^lib;
  - localization = T₄^dec − T₄^OU(x_*).
- R uses the two-model c₀ = T₄^OU(0) measured in the same runs.
- Cross-check: §11 of `decoder_moments_2026-09-29.md` found R = +0.0119, +0.0500, +0.0531 at N′ = 12, 13, 14, from
  independent realizations.

| N′ | a | samples | T₂^dec | x_* | x_exact | T₄^dec | T₄^OU(x_*) | T₄^lib | T₄^free | liberation | isotropic | localization | R_dec | R_OU(x_*) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 8 | 0.7714 | 200 | 0.8174(2) | 0.476 | 0.250 | 0.7778(2) | 0.7603(2) | 0.7627 | 0.7311 | 0.0316 (68 %) | -0.0024 (-5.1 %) | 0.0175(1) (37 %) | -0.0233 | -0.0408 |
| 9 | 0.6429 | 150 | 0.7319(3) | 0.557 | 0.333 | 0.6666(3) | 0.6373(3) | 0.6354 | 0.5609 | 0.0745 (70 %) | +0.0019 (1.8 %) | 0.0293(1) (28 %) | -0.0113 | -0.0405 |
| 10 | 0.6429 | 100 | 0.7456(3) | 0.582 | 0.400 | 0.6790(3) | 0.6469(3) | 0.6475 | 0.5609 | 0.0867 (73 %) | -0.0006 (-0.5 %) | 0.0321(1) (27 %) | -0.0174 | -0.0495 |
| 11 | 0.5260 | 60 | 0.6690(3) | 0.630 | 0.455 | 0.5799(4) | 0.5419(3) | 0.5356 | 0.4078 | 0.1278 (74 %) | +0.0064 (3.7 %) | 0.0380(2) (22 %) | +0.0138 | -0.0242 |
| 12 | 0.5260 | 40 | 0.6861(3) | 0.650 | 0.500 | 0.5964(3) | 0.5548(3) | 0.5519 | 0.4078 | 0.1441 (76 %) | +0.0029 (1.5 %) | 0.0417(2) (22 %) | +0.0121 | -0.0295 |
| 13 | 0.4248 | 16 | 0.6184(5) | 0.684 | 0.538 | 0.5101(5) | 0.4669(5) | 0.4562 | 0.2843 | 0.1719 (76 %) | +0.0107 (4.7 %) | 0.0432(2) (19 %) | +0.0515 | +0.0082 |
| 14 | 0.4248 | 14 | 0.6368(4) | 0.696 | 0.571 | 0.5283(4) | 0.4816(4) | 0.4742 | 0.2843 | 0.1900 (78 %) | +0.0074 (3.0 %) | 0.0467(2) (19 %) | +0.0508 | +0.0041 |

**3.2 Isotropic OU vs free liberation along the whole curve.**

| N′ | a | q = e^{−18/N′} | max (T₄^OU − T₄^lib) over the curve [at x] | same for T₈ | c₀ − c₀^free (two-model non-freeness) |
|---|---|---|---|---|---|
| 8 | 0.7714 | 0.105 | -0.0025 [x=0.383] | -0.0028 | -0.0010 |
| 9 | 0.6429 | 0.135 | +0.0043 [x=0.831] | +0.0099 | +0.0025 |
| 10 | 0.6429 | 0.165 | +0.0031 [x=0.000] | +0.0057 | +0.0034 |
| 11 | 0.5260 | 0.195 | +0.0081 [x=0.831] | +0.0182 | +0.0079 |
| 12 | 0.5260 | 0.223 | +0.0083 [x=0.000] | +0.0165 | +0.0080 |
| 13 | 0.4248 | 0.250 | +0.0120 [x=0.831] | +0.0258 | +0.0099 |
| 14 | 0.4248 | 0.276 | +0.0121 [x=0.000] | +0.0222 | +0.0132 |

- The largest OU − liberation deviation tracks the two-model non-freeness c₀ − c₀^free at every N′. Both grow with
  q = e^{−λ}, λ = 18/N′ (Liu–Shen: two independent models cross with weight q).
- The maximum sits at the x = 0 endpoint (even N′) or at a mid-curve bump x ≈ 0.83 (odd N′).
- In the bulk of the spectrum, OU and free liberation coincide at matched T₂ (panel d). OU has a soft tail below the
  free lower edge.



**Site-local family vs OU at equal T₂** (excess = T₄^site(z) − T₄^OU(at the same T₂); z = −1 is the decoder):

| N′ | z | T₂^site | x_*(OU) | excess T₄ | excess T₆ | excess T₈ |
|---|---|---|---|---|---|---|
| 8 | +0.924 | 0.9853 | 0.971 | +0.0001(0) | +0.0003 | +0.0006 |
| 8 | +0.707 | 0.9477 | 0.892 | +0.0021(1) | +0.0050 | +0.0082 |
| 8 | +0.383 | 0.9019 | 0.778 | +0.0068(2) | +0.0136 | +0.0191 |
| 8 | +0.000 | 0.8619 | 0.656 | +0.0121(2) | +0.0207 | +0.0261 |
| 8 | -1.000 | 0.8180 | 0.479 | +0.0174(2) | +0.0259 | +0.0303 |
| 10 | +0.924 | 0.9797 | 0.977 | +0.0003(0) | +0.0008 | +0.0014 |
| 10 | +0.707 | 0.9283 | 0.914 | +0.0039(1) | +0.0093 | +0.0153 |
| 10 | +0.383 | 0.8660 | 0.824 | +0.0121(2) | +0.0248 | +0.0358 |
| 10 | +0.000 | 0.8113 | 0.729 | +0.0213(3) | +0.0383 | +0.0507 |
| 10 | -1.000 | 0.7449 | 0.580 | +0.0324(2) | +0.0517 | +0.0637 |
| 12 | +0.924 | 0.9743 | 0.980 | +0.0005(0) | +0.0014 | +0.0027 |
| 12 | +0.707 | 0.9096 | 0.926 | +0.0057(1) | +0.0137 | +0.0221 |
| 12 | +0.383 | 0.8326 | 0.851 | +0.0166(3) | +0.0334 | +0.0478 |
| 12 | +0.000 | 0.7656 | 0.771 | +0.0279(4) | +0.0500 | +0.0663 |
| 12 | -1.000 | 0.6826 | 0.644 | +0.0414(3) | +0.0669 | +0.0835 |

**Parity family** (decoder s = 1 and parity sets |A| = s), excess over OU at equal T₂, N′ = 12 and 13:

| N′ | s | x_exact(s) | x_* | T₂ | T₄ − T₄^OU(x_*) | (T₄ − T₄^OU)/(T₄^OU − c₀) |
|---|---|---|---|---|---|---|
| 12 | 1 | +0.500 | 0.650(1) | 0.6861 | +0.0417(2) | 0.30 |
| 12 | 2 | +0.182 | 0.384(3) | 0.5766 | +0.0175(2) | 0.43 |
| 12 | 3 | +0.009 | 0.153(14) | 0.5346 | +0.0045(2) | 0.65 |
| 12 | 4 | -0.055 | -0.009(19) | 0.5252 | +0.0009(2) | nan |
| 12 | 5 | -0.045 | 0.011(18) | 0.5263 | +0.0003(2) | nan |
| 12 | 6 | +0.000 | 0.060(21) | 0.5290 | +0.0007(2) | nan |
| 13 | 1 | +0.538 | 0.684(2) | 0.6184 | +0.0432(2) | 0.25 |
| 13 | 2 | +0.231 | 0.441(5) | 0.4945 | +0.0184(2) | 0.32 |
| 13 | 3 | +0.049 | 0.249(13) | 0.4447 | +0.0061(2) | 0.37 |
| 13 | 4 | -0.035 | 0.042(25) | 0.4265 | +0.0017(2) | nan |
| 13 | 5 | -0.049 | -0.008(25) | 0.4237 | +0.0002(2) | nan |
| 13 | 6 | -0.021 | 0.044(27) | 0.4266 | +0.0004(2) | nan |


**Small-angle heterogeneity** κ = (T₄ − 2T₂ + 1)/(1 − T₂)² = E sin⁴θ / (E sin²θ)², at the point nearest the identity
(OU x = 0.981; site z = 0.924):

| N′ | free liberation 1/(1 − a) | OU | site-local |
|---|---|---|---|
| 8 | 4.38 | 4.51 (+3 %) | 5.36 (+22 %) |
| 10 | 2.80 | 3.02 (+8 %) | 3.86 (+38 %) |
| 12 | 2.11 | 2.36 (+12 %) | 3.15 (+49 %) |

**Matched OU parameter.** x_* ≈ √χ̂₂(U) = √(1 − 8(N−1)/(N(N+1))) at chaos order ≤ 1 (D5 §1). Measured x_*:
0.476, 0.557, 0.582, 0.630, 0.650, 0.684, 0.696 against predicted 0.471, 0.537, 0.588, 0.628, 0.660, 0.688, 0.711
(N′ = 8–14).
The naive per-chord x_s = 1 − 2p/N gives 0.25–0.54.

**Remarks on the parity family.**
- For multi-site probes the excess over OU is smaller in absolute terms, but larger *relative* to the OU distance from
  independence (0.30 → 0.43 → 0.65 for s = 1, 2, 3 at N′ = 12).
- At these sizes heavier parity probes are therefore not "more isotropic". This does not contradict D5's double-scaling
  statement (decoder ≡ OU for s ∝ √N with s/N → 0, p ∝ √N), which N′ ≤ 14 cannot reach.


## 4. Interpretation

**4.1 What R is.**
- The §11 split T₄ = a² + 2a(T₂ − a) + S + R isolates R as "beyond twirl data", but R is not a natural physical
  quantity. It is a small difference of two large terms of opposite sign. The natural reference is free liberation at
  the *same* T₂, which already contains the one-wormhole term and a sizeable, sign-changing remainder R_lib(η)
  (D5 §5).
- Relative to that reference, the decoder's fourth moment has two additions:
  1. an **isotropic non-freeness** (OU − liberation), small and positive at N′ ≥ 11. It grows with N′ in step with the
     two-model non-freeness S − S_free (+0.6, +2.0, +4.3 % at N′ = 10, 12, 14; progress report §5);
  2. a **localization term** (decoder − OU at equal T₂), which is the dominant non-free part.
- The localization term is a property of *which* couplings the probe touches, not of the ℤ₂ sign:
  - the Gaussian site-local family produces it continuously;
  - its ℤ₂ endpoint equals the decoder within errors;
  - resampling the site's couplings (z = 0) already gives about ⅔ of it.

**4.2 Mechanism [interp, grounded in D5 §1].**
- The site-averaged decoder weight of a chaos component is 1 − 2o(α)/N, where o is the odd cover of the cross-chord
  index sets. Isotropic OU weighs the same component x^{|α|}. The two agree exactly when cross-chord index sets are
  independent.
- In N=2 SYK, {Ψ_J, Ψ_K†} = 0 for disjoint J, K (odd p). So the non-trivial chaos content of P is built from
  *overlapping* index-set pairs: e.g. P₍₁₎ ∝ H^(2) (D4), whose terms have |J∩K| ≥ 1. Chords through site i therefore
  cluster. The number of probe-sensitive chords is then overdispersed relative to an isotropic deformation. At fixed
  mean (equal T₂) this raises the variance of sin²θ (T₄ − 2T₂ + 1 = E sin⁴θ).
- Spectrally (panel d), it is the wall pile-up at cos²θ ∈ {0, 1}: BPS states in which site i is effectively frozen or
  maximally uncertain. This is the UV imprint of n_i on the LMRS operator Pn_iP (cf. CCMP Fig. 5a), absent for
  isotropic deformations.
- Infinitesimally it is the excess normalized variance of the site-restricted quantum-metric operator (κ_site vs
  κ_OU ≈ κ_free = 1/(1 − a)).

**4.3 Gravitational reading [interp].**
- **Isotropic part.**
  - The OU four-arc function is, in double scaling, the overlap of two mismatched two-sided BPS wormholes (left/right
    boundaries with couplings ω, ω_x; D5 §3), i.e. an out-of-time-order correlator of the q-OU process.
  - At N′ ≤ 14 (q = e^{−18/N′} ≈ 0.1–0.3) it is close to free liberation. That suggests, as a conjecture, that the
    double-scaled isotropic function is a *q-deformed liberation* that reduces to the free one as q → 0, consistent
    with Liu–Shen's q₁₂ = q.
  - This part is universal and is what a two-colour chord computation in BBNR/BLY language should reproduce.
- **Localization part.**
  - In strict double scaling a single-site probe trivializes (x_s → 1). For heavy probes (s ∝ √N) the index-set
    correlations are subleading, so decoder ≡ OU (D5 §1). The localization term is therefore a finite-p, site-level
    (UV) effect.
  - It is *not* captured by structureless matter chords (BINT/BLY random matter) or by super-JT. Any gravity
    computation of the decoder's R must either include site-resolved chord rules or explicitly exclude this term.
- **Consequence for action item 2 of `progress_2026-09-29.pdf` §7.3.**
  - The four-arc chord computation should target the *isotropic* function T₄^OU(x), which is free-liberation-like and
    universal, plus a separate finite-p correction.
  - Expecting a single chord calculation with a "ℤ₂ flag per chord" to reproduce the decoder's R at p = 3 is not
    justified by these data.

## 5. Status and open items
- **Exact [derived]:** D5 §§1–5 (generating function; four-arc rule; Krawtchouk weights; OU sum rule; small-angle
  QGT; conjugation flow; liberation closed form).
- **Numerical [num]:**
  - OU ≈ free liberation at equal T₂ (deviations ≲ 0.012 in T₄ at N′ ≤ 13, growing with N′);
  - decoder localization excess, its site-family origin, and its infinitesimal (κ) signature;
  - the three-way split.
- **Interpretation [interp]:** mechanism (chord clustering), q-deformed liberation, the gravity division of labour.
- **Open.**
  1. An analytic model of the localization term from the site-restricted quantum geometric tensor.
  2. N′ = 15–16 for the trend of its share.
  3. The two-colour chord computation of T₄^OU in double scaling (BBNR/BLY), and a test of "q-deformed liberation".
  4. The Krawtchouk (u_t) grading of the parity family.
  5. Whether the isotropic non-freeness has a closed form in terms of S − S_free.

## 6. Reproducibility
- Commit: see git log (branch `t4-fourarc`). Runs used code at `1aa3dd3-dirty` (recorded in each `meta_*.json`).
- Commands (all under `scripts/memwatch.py`):
  - OU: `run_ou_fourarc.py --Np N --samples S --seed0 0` for (N, S) = (8, 200), (9, 150), (10, 100), (11, 60),
    (12, 40), (13, 16); K = 16, probe-proj 2, parity-sets 3.
  - N′ = 14: `--samples 14 --fast-basis --pair-starts 4`.
    - The fast basis comes from Gram-matrix eigh + QR (BLAS3). It spans the same subspace as the Hodge/SVD basis
      (min overlap singular value 1 − 3·10⁻¹⁵ at N′ = 10, 12) and reproduces a reference N′ = 10 sample exactly.
    - With `--pair-starts 4` only the pairs starting at i ∈ {0, 4, 8, 12} are used, 4 per lag instead of 16; pairs at
      the same lag within a sample are strongly correlated anyway.
    - A first N′ = 14 attempt with the SVD basis (started 18:12:24) was stopped after about 5 min, before completing a
      sample. A stack sample showed it in single-threaded `zgesdd`/`zgemv`. Its meta file is kept as
      `aborted_meta_p3_N14_F7_20261005T181224.json`.
    - The replacement uses the Gram-eigh basis (5.8 s per projector at N′ = 14, measured) and `eigvalsh(OO†)` pair
      spectra (1.1 s vs 2.3 s for `svdvals`; equal to 10⁻¹⁴).
    - The old basis was not benchmarked at N′ = 14, so no speed-up factor is claimed.
  - Site family: `--family site --seed0 1000 --probe-proj 1` for (8, 150), (10, 80), (12, 30).
- Logs: `log_small.txt`, `log_large.txt`, `log_site.txt`, `memwatch_*.jsonl`.
- Runtimes (three jobs concurrently on 10 cores): N′ = 12 ≈ 20 s, 13 ≈ 70 s, 14 ≈ 212 s per sample (fast basis, run alone).
- Peak memory: 0.37 GB at N′ = 13 and 1.43 GB at N′ = 14 (memwatch); the figure script peaks at 0.44 GB.
- Analysis: `analyze_ou_fourarc.py --Np 8 9 10 11 12 13 14 --site --markdown --json …/analysis.json`.
- Figure: `plot_ou_fourarc.py --Np 8 9 10 11 12 13 14 --panel-a 12 14`.
