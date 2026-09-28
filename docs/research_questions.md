# Research questions

*Living document. Last revised 2026-09-28 after the project archaeology audit (`docs/project_archaeology.md`).*
Status tags: **[open]**, **[in progress]**, **[answered]**, **[superseded]**.

## Central question

> The uplift decoder O_T = P_B Π_T P_B of one-flavor N=2 SYK (Π_T = 1 − n_{N′}, the "slot") has a GUE-correlated
> spectrum on [0,1] (Chryssanthacopoulos–Vegh 2026). **What law governs its spectral density at finite and large N,
> and can that law be derived from N=2 super-JT gravity (LMRS-type projected-operator calculus) rather than
> postulated?**

**Status 2026-09-28 (`docs/derivations.md` D2).** The second moment is now derived gravitationally. The
finite-λ N=2 super-chord theory (Boruch–Lin–Yan) gives r with no free parameter, to 2 %. Structurally,
r = a + Σ_{n≥1} P_n q^{2Δn}: the Wachter/Haar value a is the zero-length-wormhole term, and the deviation from
freeness is the finite-length-wormhole contribution. Predicted: η_r peaks near N′ ≈ 22. Open: higher moments / the
density (multi-particle chord wormholes), X and e₃.

The previous draft answered "Wachter/MANOVA (two free projectors) + a Haar genus expansion = gravity". The audit
found that this answer describes a **Haar null model**, not the SYK decoder: its gravitational derivation is not
valid, and the SYK spectrum departs from it at O(1). The questions below replace that framing.

Notation: N = old size, N′ = N+1, P = ⌊N′/2⌋ (down channel), D = C(N′,P), d = dim B^P_{N′}, m = C(N,P),
a = d/D, b = m/D. B⁰ := B^P_N ⊕ (B^{P−1}_N ∧ e_{N′}) is the LES reference subspace.

---

## Q1. Is the decoder variance exponentially small (Haar) or power-law small (LMRS)? **[in progress — Haar excluded; finite-λ chord prediction matches to 2 %; asymptotics open]**

**Update 2026-09-28 (`docs/derivations.md` D1).** The LMRS zero-energy two-point function in our single R-charge
sector (j = 0 for even N′, −1/6 for odd N′) gives a parameter-free prediction r_LMRS ∝ N′^{−2/3}.
- As derived, it **overshoots** the measured q=3 r by 15–24 % (N′ = 9…15), with the gap shrinking in N′.
- Its j-factor alone correctly accounts for the even/odd alternation (parameter-free).
- A fitted correction, r ≈ r_LMRS (1 − 2.2/N′), fits well, but a plain power law N′^{−0.44} fits equally well.
  The pure N′^{−2/3} shape with free normalization does not fit (χ²/dof 335).
- LMRS themselves compare only at N = 16, with no finite-N correction, and are ~10 % off for their neutral bilinear.
- All Haar-based models fail.

**Out-of-sample R-charge test (D1, pre-registered).** At fixed even N′, the j = ±1/3 sectors are predicted to have
Var_{1/3}/Var_0 = 0.645 (no free parameter). Measured: 0.642, 0.668, 0.677, 0.683 (N′ = 8…14). The Haar null gives
0.56–0.59, and no j-dependence gives 1. So the prediction holds to within 6 %, but the data drift away from it with
N′, not toward it.
**Resolution at finite λ (D2, 2026-09-28).** The double-scaled super-chord two-point function of Boruch–Lin–Yan
(eq. 4.8), evaluated at λ = 2p²/N′ = 18/N′ and Δ = 1/3 with no free parameters, reproduces the measured r to within
2.2 % for N′ = 5…16 (both parities). The super-Schwarzian overshoots by 15–60 %. So at our sizes the decoder variance
is in the double-scaled regime (λ ≈ 1.1–3.6), not the Schwarzian one. Residual: a drift of −2 % → +2 % (j = 0,
∓1/6) and up to +7 % in the j = ±1/3 sectors at N′ = 14. Since DSSYK's λ → 0 limit is not the fixed-p = 3 large-N
theory (α_S = 1/72 vs 0.00842), the agreement must fail eventually. Q1 status: Haar excluded; the variance is
quantitatively a chord/gravity quantity at accessible N′. The asymptotic N′ → ∞ power remains open.
**Why they differ (D1, "no conformal window").** The zero-energy correlator samples boundary times
Jt ~ 2α_S N′ ≈ 0.75 at N′ = 15, below the microscopic time 1/J. Half its weight comes from where the conformal input
exceeds the true correlator's UV maximum, so the Schwarzian formula is used outside its regime, and the overshoot is
expected. A no-parameter "saturated conformal" toy reproduces the size and trend. The N′^{−2/3} regime needs
N′ ~ 10² (median Jt ≈ 10 at N′ ≈ 200).
Open: insert the exact large-N Schwinger–Dyson correlator into the zero-energy average (controlled finite-N′
prediction); the same issue applies to predictions of X and e₃.

**Measured 2026-09-28** (`scripts/audit/variance_scaling.py`, seeds 0–3, 2 for N′=14; half filling P=⌊N′/2⌋):
r ≡ Var(λ)/[b(1−b)] (fraction of Π_T's UV variance surviving BPS projection; Haar/Wachter predicts r = a):

| N′ | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
|---|---|---|---|---|---|---|---|---|
| a | .771 | .771 | .643 | .643 | .526 | .526 | .425 | .425 |
| r | .797 | .818 | .735 | .736 | .667 | .699 | .612 | .652 |
| r/a | 1.03 | 1.06 | 1.14 | 1.14 | 1.27 | 1.33 | 1.44 | 1.53 |

The deviation from Wachter **grows** with N (δm₂ = 0.006 → 0.057), i.e. it is not a finite-size correction to
Wachter. r decays far more slowly than a. Power law vs slow exponential cannot be decided from 8 points.

**Chord crossing weight (2026-09-28, `scripts/audit/crossing_ratio.py`).** For two projected occupations
Â = P_B(n−⟨n⟩)P_B of different modes, X = τ(ÂᵢÂⱼÂᵢÂⱼ)/τ(ÂᵢÂᵢÂⱼÂⱼ). X→0 is free (non-crossing, Wachter),
X→1 is a light/commuting probe.

| N′ | 8 | 9 | 10 | 11 | 12 | 13 |
|---|---|---|---|---|---|---|
| X_SYK | .851 | .800 | .801 | .766 | .792 | .735 |
| X_Haar | .777 | .651 | .652 | .533 | .531 | .431 |

X_Haar ≈ a, whereas SYK crossings stay O(1) and drift slowly. This is the light-probe chord regime, and X is a
directly gravity-computable number (disk 4-point OTOC/TOC at Δ = 1/3).

**Chord-invariant campaign (2026-09-28; `research/notes/chord_invariants_2026-09-28.md` §3).** Exact, q=3
N′ = 5…15 and q=5 N′ = 9…16, with a Haar null model that reproduces free compression. Findings:
(i) r/a grows to 1.70 at N′=15, with r ≈ 0.89 a^0.41 (asymptotic form undetermined).
(ii) Low free cumulants look like Wachter(t ≈ r), but the spectra are not Wachter(r) (KS 0.15–0.23).
(iii) F = 1, and words with ≤ 2 crossings are multiplicative in X.
(iv) The fully crossing word scales as X^{e₃} with **e₃ = 2.27 → 2.41 (± 0.01)**, drifting upward: neither free
(2) nor the single-parameter chord rule (3).
Point (iv) is the main quantitative target for a gravity computation.

**Earlier estimates.** SYK variance exceeds Wachter by 0.02–0.035 at N′ = 9–11 (audit smoke test) and δm₂ ≈ 0.035–0.044 at
N = 10–12 (draft, unverified code). A crude, *unverified-normalization* evaluation of LMRS eq. (66) (j = 0,
Δ = 1/3, C = 0.00842 N) gives O(0.1) at N′ = 13.

**To do.**
1. ~~Analytic: derive the zero-energy two-point function of O = n_{N′} − ⟨n⟩ in the fixed-p BPS sector.~~ Done
   2026-09-28, `docs/derivations.md` D1 (i = k handled at leading order; j per parity; b_q̂, α_S from LMRS).
   Issues: i = k bilinear vs LMRS's i ≠ k; R-charge sector j for our P; odd N′ shift (LMRS fn. 6); normalization
   b_q̂, α_S. Output: prediction V_LMRS(N).
2. Numerical: Var(λ) and m₃ − … (connected moments) vs N for N′ = 8…16 (exact) and larger with `src/q_scan_mf.py`
   (Hutchinson). Plot against ab(1−b) and against c·N^{−2/3}.
3. Compare also with the "UV imprint" picture (Berry Fig. 5a).

**Success criterion.** A fit that distinguishes exponential from power-law decay with error bars, and an analytic
prediction matching within LMRS-level accuracy (LMRS Table 1 achieved ~10%).

## Q2. Is the decoder a partially liberated version of the LES reference subspace? **[open — promising]**

**Motivation.** CV eq. (27) implies d = dim B^P_N + dim B^{P−1}_N (verified in all repo data, N = 8…16). So B⁰ has
the same dimension as the enlarged BPS space B, and O_T restricted to B⁰ is exactly Bernoulli on {0,1}.
Exploratory check (`scripts/audit/reference_subspace_check.py`, N′ = 10–12, 2 seeds each):
- mean cos² of principal angles between B and B⁰ ≈ 0.81–0.88; 55–65% of angles have cos² > 0.9;
- eigenvectors with λ < 0.05 carry 97–99.6% of their weight in B^{P−1}_N ∧ e; those with λ > 0.95 carry ~99% in B^P_N;
- bulk eigenvectors (0.2 < λ < 0.8) carry only ~50–58% in B⁰.

**Hypothesis.** B = rotation of B⁰ generated by the new couplings. The decoder spectrum interpolates between
Bern(dim B^P_N/d) (no rotation) and Wachter(a,b) (full liberation, cf. free Jacobi process, Demni–Hamdi–Hmidi).
The wall pile-up and "soft tails" are the un-liberated remainder. This replaces the unsupported "near-BPS energy
image" explanation.

**To do.** Perturbation theory in the new couplings D (B = B⁰ + O(D) corrections via the LES / Hodge Green
function); the distribution of principal angles between B and B⁰ vs N; test whether a free-Jacobi-process law at
some effective time t(N) fits the SYK spectrum; check whether t(N) grows (→ Wachter) or saturates.

## Q3. Exact endpoint atoms: generic-position counting **[partly answered]**

**Known.** n₁ = dim(B^P_N ∩ ker D∧), n₀ = dim(B^{P−1}_N ∩ ker ι_D̄) (exact). Fractions: 55.6% (N = 8, 9), 32.9%
(10, 11), 5.6% (12, 13), 0 (14, 15, 16). The counts nearly saturate the rank lower bounds
n₁ ≥ dim B^P_N − C(N,P+2), n₀ ≥ dim B^{P−1}_N − C(N,P−3) (exact at N = 8, 9; small excess at N = 10–12; +27 at N = 13).
Since dim B ~ 3^{N/2} and binomials ~ 2^N, the atoms must vanish at large N. This is **not** a resonance.

**To do.** Explain the excess by generic position inside CV-type constrained rooms (CV Table 1 logic). Prove the
vanishing for N ≥ 14 conditional on generic position. Update the draft's narrative.

## Q4. What is a legitimate gravitational description of O_T? **[open]**

**Requirements from the literature.** (i) e^{S0} ∝ BPS degeneracy d, not D (LMRS eq. 82). (ii) Π_T enters as a
boundary operator insertion (light neutral bilinear, Δ = 1/3), not as an EOW-brane species (no derivation exists).
(iii) BPS projection does *not* trivialize the super-Schwarzian: zero-energy correlators depend on Δ and R-charge.
(iv) Disk crossings are weighted by Δ-dependent OTOC factors, O(1) in e^{S0}.
(v) Boruch–Iliesiu–Yan (2023) is the worked BPS analogue of PSSY. Matter-prepared BPS states have Marchenko–Pastur
Gram statistics with dimension e^{S0}cos(πj). Trumpet/higher-genus contributions **vanish** at β→∞, and the
statistics are exactly those of Haar-random vectors *inside* H_BPS. Gravity therefore predicts ETH-like random
matrix elements ⟨ψ_i|Π_T|ψ_j⟩ with LMRS variance, not Haar randomness of the BPS subspace in the Fock space and
not a 1/D² handle tower. (vi) Finite-Δ crossing weights from chords: q̃ = q^Δ (Berkooz–Mamroud eq. 2.28); the
supersymmetric chord methods of Boruch–Lin–Yan (`papers/boruch_lin_yan_2023.pdf`, noted 2026-09-28) are the
computational route. Their §5 builds one-particle wormholes and the zero-T OTOC operator (eq. 5.18) without
evaluating it. With site-local crossing rules, that is a direct finite-λ prediction of X and e₃ (see D2).

**To do.** Compute m₂ (Q1), then the disk 4-point function at finite Δ (LMRS eq. 145-type integral) to get m₃, m₄.
Double-scaled SYK chord rules may give finite-Δ crossing weights. Decide whether any E=0 topological model with
e^{S0χ} weights is appropriate (it would describe the Haar null model only).

## Q5. The soft tail: correct formulation and large-N fate **[open]**

**Known.** A reproducible O(10%) population below λ₋ at N = 12–15 (1–3 seeds). Haar predicts exponentially few
(Collins). The exact identities λ = E_b/(E_b+σ_a²) and 1−λ = E_a/(E_a+σ_b²) hold but are near-tautological.
Tail vectors are ≈ 99.6–99.8% harmonic in their occupied component (prior audit).

**To do.** (a) Derive and numerically verify the correct constrained generalized eigenproblem for O_T in terms of
old-theory data (the prior audit's proposal). (b) Define a threshold-independent tail measure (the threshold
λ₋(a) moves with N), e.g. spectral weight within distance δ of the walls, or the B⁰-overlap distribution of Q2.
(c) Larger N (N′ = 16–18) with sparse methods. (d) Compare edge fluctuations with the Tracy–Widom scale of the
matching Jacobi ensemble.

## Q6. BPS chaos of the decoder: strength, Thouless time, and relation to fragility **[open]**

**Known.** ⟨r⟩ = 0.599 ± 0.001 at N′ = 12, P = 6 (CV). Tail ⟨r⟩ ≈ GUE (draft, provisional, no code).

**To do.** SFF of O_T in **SYK** (not a Jacobi surrogate), with Gaussian filtering. Thouless time vs N. Test CV's
conjecture "metric fragility ⇔ BPS chaos" on (i) Miyahara–Shibuya's integrable interpolation Q_g (fortuitous but
integrable), (ii) the two-flavor model (monotone vs fortuitous), (iii) the single-matrix model (rigid).

## Q7. Which decoder spectrum is the physical object? **[open, clarification]**

CV's diagnostics are state-dependent measures dμ_{α,T} for *old* BPS states α. The draft uses the uniform
spectrum of O_T on the *enlarged* BPS space. `src/KT_histogram.py` computes a third object, the compression of
K_T onto B^p_N (eigenvalues of B_old† P_{B_new} B_old). Decide which object each claim refers to, and derive the
free-probability null model for each (the Wachter overlay in KT_histogram.py is not derived).

## Q8. Dependence on q and on double scaling **[open, low priority]**

`src/q_scan.py` / `q_scan_mf.py` test whether δm₂ falls with q at fixed (a,b). LMRS intuition suggests the
opposite for a *bilinear* probe (Δ = 1/q̂ gets lighter). Worth running once Q1 fixes what "freeness" should mean.

## Q9. Entanglement reading **[open, low priority]**

H(λᵢ) is exactly the entanglement entropy of decoder eigenvector i across the single-mode cut. Study
(1/d)Σ H(λᵢ) vs the Wachter value and vs the Haar-typical-state value H(b). The draft's "fragility = entanglement
deficit" is a restatement of the wall pile-up and should be framed that way.

---

### Superseded questions (kept for the record)

- "Does the super-JT action reproduce the Weingarten handle tower to all orders?" — **[superseded]**: the handle
  tower belongs to the Haar null model. The gravitational identification failed (e^{S0} = D, EOW-brane slot,
  trivialized Schwarzian).
- "Is the soft tail the image of the smaller theory's near-BPS density ρ_H?" — **[superseded]** by Q2/Q5: tail
  eigenvectors are near-harmonic, not near-BPS-energy eigenstates.
- "Is the protected-atom sector a finite-N resonance at N ≈ 12–13?" — **[answered, no]**: see Q3.
