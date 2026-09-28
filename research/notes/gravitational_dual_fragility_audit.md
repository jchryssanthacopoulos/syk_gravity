# Claim-by-claim audit of the current draft

**Document audited:** `research/tex/gravitational_dual_fragility.tex` (= `research/pdfs/gravitational_dual_fragility.pdf`),
"Fragile Fortuity: The Gravity of BPS Overlaps in Supersymmetric SYK".
**Prior audit:** `research/pdfs/decoder_notes_audit.pdf` (Claude-generated, Sept 2026). Where this note agrees
with it, the point was re-checked independently. Where it disagrees, that is stated explicitly.
**Date:** 2026-09-28.

Legend: ✅ verified · ⚠️ correct but mislabelled or overstated · ❌ incorrect or unsupported · ❓ not checkable from repo.

## 0. What the draft is trying to say (in one paragraph)

The draft takes CV's decoder O_T = P_B Π_T P_B (one-flavor N=2 SYK, uplift N → N+1, slot = "new mode empty")
and argues three things. (i) Its disorder-averaged spectrum is the Wachter law of two free projectors. (ii) The
1/D² genus expansion of that law "is" a sum over topologies in the E=0 sector of N=2 super-JT with two species
of EOW branes. (iii) What escapes this description (exact λ=0,1 atoms and a ~10% "soft tail" below the Wachter
edge) is the image of the smaller theory's near-BPS spectrum. Around this it builds chord-diagram, topological
recursion, Kunisky-error and bordered-Toeplitz machinery.

## 1. Section 1 — decoder side

| # | Claim (draft location) | Status | Notes |
|---|---|---|---|
| 1.1 | O_T = P_B Π_T P_B, 0 ⪯ O_T ⪯ 1, λ = cos² of principal angles (eq. decoderdef) | ✅ | Same as CV eqs. 101–102. |
| 1.2 | F_T = (1/d)Tr O_T = m₁, κ_T = m₋₁ − 1 | ⚠️ | **Redefinition.** In CV, F_T(α), κ_T(α) are state-dependent for *old* states α (CV eqs. 50, 61). Also m₋₁ = ∞ when exact λ=0 atoms exist (audit); finite values need a pseudo-inverse or restriction. |
| 1.3 | Hodge projector P_B = 1 − QGQ† − Q†GQ (eq. PB) | ✅ | Standard. |
| 1.4 | Four block BPS conditions (eq. bpsfour) | ✅ (up to signs) | Consistent with CV eqs. 20–22, 36–39, 96. |
| 1.5 | b = −G Q_N† (D∧a) (eq. belim) | ❌ | Omits the harmonic part: general solution b = b₀ − GQ†u with b₀ ∈ ker H_N. Exactly the λ=0 atoms (b ∈ B^{p−1}_N) are excluded. |
| 1.6 | O_T = (1 + C′†GC′)^{−1} "on ker Q_N" (eq. KTself) | ❌ | (a) harmonic b₀ omitted; (b) domain wrong: a must also satisfy D∧a ∈ im Q_N and condition (iv); (c) even on the right domain it is a *generalized* eigenproblem, [P_A(1+Σ)P_A]^{−1}, not a compressed inverse. The formula cannot produce λ = 0. Audit gives a constrained generalized-eigenproblem replacement (unverified here). |
| 1.7 | Feynman-diagram / Dyson-series reading (Sec. 1.3, Fig. 1) | ❌ as derived | Rests on 1.6. |
| 1.8 | Disorder-averaged spectrum is Wachter "because C′ is independent Gaussian and B is (numerically) Haar-typical" | ⚠️ | Assumption, not derivation. B depends on both old and new couplings. Figures in `results/figures/decoder_wachter*.png` show a visibly non-Wachter interior (mid-band depleted by ~2×, walls enhanced); KS 0.118 (12→13), 0.167 (14→15). |
| 1.9 | Wachter density, edges, atoms, S-transform, moments m₁–m₃, spectral curve (eqs. wachter–curve) | ✅ | Checked by hand and in `scripts/audit/weingarten_W1_check.py` (disk terms). |
| 1.10 | "⟨λ⟩ = b is an exact trace identity at every N" | ⚠️ | True for the **disorder average** (exchange symmetry of modes + Σᵢnᵢ = p). Not per realization: measured m₁ = 0.536–0.560 at N′=9 (b = 0.556). Per realization only at even N′ and half filling (particle–hole ⇒ m₁ = 1/2 exactly). |
| 1.11 | KS ≲ 0.1 "tightening with N" (N = 6…10) | ⚠️ | Audit reports non-monotone KS 0.111, 0.060, 0.112, 0.046, 0.101 (even/odd N′ alternation likely). |
| 1.12 | Protected atoms: λ=1 ⇔ a ∈ B^p_N ∩ ker(D∧); λ=0 ⇔ b ∈ B^{p−1}_N ∩ ker ι_D̄ | ✅ | Matches CV/audit and `src/protected_atoms.py`. |
| 1.13 | Atoms are a "finite-N resonance … appreciable only where dim B^p_N ≈ C(N,p+2), near N=12–13" | ❌ | Rerunning `protected_atoms.py` (N = 8…11, seed 0) gives atom fractions **55.6%, 55.6%, 32.9%, 32.9%**, then 5.6% (12, 13) and 0 (14–16). This is a monotone decay, not a resonance. The counts nearly saturate the rank bounds n₁ ≥ dim B^p_N − C(N,p+2) and n₀ ≥ dim B^{p−1}_N − C(N,p−3) (exact at N = 8, 9; excess +2 at 10–12; +27 at 13). |
| 1.14 | d = 3^{⌊N/2⌋} (×2 for odd N), "two-flavor tower" dim B^m_{2m+1} = 3^m "proved by generating-function/Lagrange–Bürmann collapse" (eq. index) | ⚠️ | Dimensions correct (N = old size; target N′ = N+1). No proof is present in the repo. CV App. A proves the walk count conditional on Conj. A.1. "Two-flavor" is the wrong model name; S₀ here means BPS entropy, clashing with S₀ = ln D elsewhere. |
| 1.15 | Split fragility = Wachter + atom resonance (eq. split) | ❌ | Omits the macroscopic non-Wachter redistribution (Sec. 4 of the draft itself). |

## 2. Section 2 — genus expansion (Haar null model)

| # | Claim | Status | Notes |
|---|---|---|---|
| 2.1 | Weingarten formula (eq. weingarten) | ✅ | Equivalent to Σ Wg(στ⁻¹) m^{c(γσ)} d^{c(τ)}. |
| 2.2 | m_k^{(1)} for k = 2, 3, 4 (eqs. m2handle–m4handle); m₁ has no correction | ✅ | Verified exactly for k ≤ 5 (this audit) and k ≤ 6 (prior audit). |
| 2.3 | Numerical check "at a=b=1/2, m₂ = 3/8 − ¼D⁻²" | ❌ | Correct: 3/8 − ⅛D⁻² (= −b(1−a)(1−b)). m₃, m₄ values quoted are correct. |
| 2.4 | Exact m₂ = b[D²(a+b−ab) − 1]/(D²−1) (eq. m2exact) | ✅ | |
| 2.5 | Non-crossing partitions / Narayana / free cumulants of Bern(a) (eq. ncmoments) | ✅ | Standard free probability. The "chord encloses a run at no cost" wording is heuristic. |
| 2.6 | "Non-crossing limit in LMRS/Berry is Δ→∞; our Δ-role is played by D→∞" | ❌ | LMRS/Berry: disk crossings are controlled by Δ and are O(1) in e^{S0} at finite Δ. Haar freeness is a different mechanism. |
| 2.7 | "One crossing = one handle ∝ 1/D²" (Fig. cross) | ❌ for gravity, ⚠️ for Haar | True only as bookkeeping for the Haar Weingarten sum. |
| 2.8 | Cylinder W_{0,2} (eq. cylinder) | ✅ (standard form) | Universal one-cut two-point function; normalization (1/d vs 1/D) not addressed. |
| 2.9 | SFF ramp t/2π, plateau d at t_H = 2πd (eq. ramp, Fig. sff) | ⚠️ | Checked only on a **Jacobi ensemble** (a=b=½), not SYK. Slope t/2π coincides with band width 1 at a=b=½. In general, SFF_c ≈ (t/2π)×(band width) before the plateau. No code in repo. |

## 3. Section 3 — "gravitational model"

| # | Claim | Status | Notes |
|---|---|---|---|
| 3.1 | Two EOW species: ψ-branes (d BPS flavours), t-branes (m slot flavours) | ❌ unsupported | Slot states are non-BPS Fock states. No boundary condition for Π_T is derived. |
| 3.2 | "BPS states at E=0 ⇒ super-Schwarzian frozen onto δ(E): only topology" | ❌ | Contradicts LMRS (zero-energy correlators depend on Δ, j). |
| 3.3 | Amplitude table: disk = e^{S0} = D, wormhole 1/(D²−1), loops d, m | ❌ as gravity / ⚠️ as Haar | e^{S0} ∝ BPS count d (LMRS eq. 82), not D. 1/(D²−1) is Wg(e) at k=2 only. |
| 3.4 | 't Hooft mechanism: planar wormholes unsuppressed, Schwinger–Dyson = subordination | ⚠️ | Correct statement about the Haar/Weingarten sum; gravitational reading is interpretation. |
| 3.5 | MP law is the b→0 limit | ⚠️ | Requires rescaling λ/b (audit). |
| 3.6 | "Four bulk surfaces of m₂" with weights b², −ab², ab, −b/D² | ✅ as Weingarten terms | Negative weights have no geometric interpretation. Call them "Haar fatgraphs". |
| 3.7 | m₃: "12 planar, 21 genus-one, 3 genus-two surfaces" | ❓ | Not re-checked; audit says valid as permutation-term classification. |
| 3.8 | W₁(x) = −b(1−a)(1−b)x(x−1)/[(x−λ₋)(x−λ₊)]^{5/2} (eq. W1) | ✅ (Haar model) | Verified through k = 5 (here) and 6 (audit), including m₅^{(1)}. "Resums the whole handle tower" is wrong: it is genus one only. |
| 3.9 | Normalization "−1/a², fixed at six (a,b) points" | ⚠️ | Should be derived. Plausibly (D/d)² from expanding in 1/d² vs 1/D², up to sign convention. |
| 3.10 | x(x−1) zeros ⇒ "handle never touches the protected atoms" | ⚠️ | Algebraic fact about the Haar W₁. The Haar model has only forced atoms; physical excess atoms are outside it anyway. |

## 4. Section 4 — "amplitude dictionary from the super-JT action"

| # | Claim | Status | Notes |
|---|---|---|---|
| 4.1 | Topological term gives e^{S0χ} | ✅ trivially | |
| 4.2 | "Decoder is built entirely from zero-energy states ⇒ dressing trivializes" | ❌ | Π_T is not built from zero-energy states. LMRS show projected correlators are dynamical. |
| 4.3 | Fatgraph duality ⇒ gravity sum "identically" = Weingarten average | ❌ | Haar genus expansion involves Weingarten/monotone-Hurwitz structure with signs, not a plain D^χ fatgraph sum. No gravitational derivation is given. |
| 4.4 | "Derived (exact) … all-genus reproduction" (Status paragraph) | ❌ | Should be "conjectural interpretation of the Haar null model". |

## 5. Section 5 — soft tail

| # | Claim | Status | Notes |
|---|---|---|---|
| 5.1 | Sub-edge population exists, O(10%), non-Haar | ✅ numerically | `results/data/softtail_scaling.csv`: N=14: 15.3–15.9%; N=15: 13.7–14.9%. N=12, 13 values (9%) are not in repo data but reproduced by audit at N′=13 (9.5–10.5%). Haar: exponentially rare (Collins). |
| 5.2 | "Reproduces every positive moment m_{k≥1}" (opening of Sec. 4 of the draft) | ❌ | Impossible on [0,1] with a macroscopic tail, and contradicted by the draft's own Table (m₂ = 0.431 vs 0.396). |
| 5.3 | Tail "grows" with N (9% → 15%) | ⚠️ | The threshold λ₋ moves with a(N); 4 sizes, 1–3 seeds. Growth vs saturation undetermined. |
| 5.4 | Bulk KS 0.08 → 0.07 "universal and improving" | ⚠️ | Conditional on restricting to [λ₋, λ₊]. The full interior histogram is visibly non-Wachter. |
| 5.5 | κ values 115, 282, 44, 45 | ⚠️ | Infinite if atoms included. Here computed on λ ≥ 1e−9 (code). |

## 6. Section 6 — "gravitational origin of the soft tail"

| # | Claim | Status | Notes |
|---|---|---|---|
| 6.1 | λ = E_b/(E_b + σ_a²) (eq. identity) | ✅ exact | But near-tautological: E_b and σ_a² are read from the same eigenvector. Machine-precision agreement is automatic. |
| 6.2 | Mirror 1 − λ = E_a/(E_a + σ_b²) (eq. identityup) | ✅ exact | Same caveat. |
| 6.3 | "The linear law is a theorem" | ❌ | The identity is a theorem. Constancy of σ_a² (hence λ ≈ cE) is an empirical correlation (0.93–0.95). |
| 6.4 | ρ_tail(λ) = c⁻¹ρ_H(λ/c): tail = image of near-BPS spectrum (eq. tailimage) | ❌ | E_b is a Rayleigh quotient, not an eigenvalue. Audit: tail vectors have 99.57–99.80% of ‖b‖² in B^{p−1}_N, with Ē⊥ ~ 16–29 (ordinary). This audit's check (`scripts/audit/reference_subspace_check.py`, N′ = 10–12): eigenvectors with λ < 0.05 have 97–99.6% weight in the occupied copy of B^{p−1}_N; λ > 0.95 have ~99% in the empty copy of B^p_N. The tail is a **near-intersection** phenomenon with the LES reference subspace, not a near-extremal-throat image. |
| 6.5 | Figures sourcing.pdf etc. | ❓ | No code or data in repo. |

## 7. Section 7 — Kunisky error terms

| # | Claim | Status | Notes |
|---|---|---|---|
| 7.1 | Citation Kunisky Lemma 4.6, Def. 4.2 | ✅ | Accurate. |
| 7.2 | Tr(O^k) = D m_k^W + Δ_k (eq. errdef) | ❌ typo | Should be d·m_k^W (m_k^W per-BPS normalized). The literal version is off by (D−d)m_k^W (checked numerically). |
| 7.3 | K_{k,ℓ} closed form (eq. Kclosed) via CS decomposition / Dickson polynomials | ✅ | Identity verified numerically for k ≤ 7 on Haar, coordinate and strongly aligned projector pairs (~1e−14). |
| 7.4 | "Verification … machine precision" (table) | ⚠️ | **Tautological**: the identity holds for *any* pair of projections of those ranks. No physics is tested. |
| 7.5 | Chebyshev / hopping-chain / "bound states" interpretation | ⚠️ | Correct mathematics, interpretive. |
| 7.6 | Decoder is "length-two bordered Toeplitz"; closed-form ρ_dec with quartic P₄ (eqs. borderedtoeplitz–P4) | ❌ / ❓ | A finite-length bordered-Toeplitz law has a.c. support exactly [λ₋, λ₊] plus finitely many atoms. It **cannot** represent a continuous sub-edge tail, contradicting Sec. 5. "Toeplitz after one step" and P₄ unverified (no code). |

## 8. Section 8 — large N

| # | Claim | Status | Notes |
|---|---|---|---|
| 8.1 | d(N) = 3^{N/2} (N even), 2·3^{(N−1)/2} (N odd), N = old size | ✅ | Consistent with CV Table 1, Berry Table 2, repo data. **The prior audit's claim that these are "interchanged by 3/2" is itself a labelling error** (it uses target N). |
| 8.2 | a ~ ½√(πN/2)(√3/2)^N | ⚠️ | Correct for even N; odd N has an extra factor 2/√3. |
| 8.3 | Wachter band collapses to δ(λ − ½), width ≈ 2√a | ✅ (Haar model) | |
| 8.4 | "The bulk disappears; physics migrates to what the sum misses" | ❌ | In the Haar model the spike carries all weight (audit). Real question: is the SYK variance exponentially or power-law small? (research Q1) |
| 8.5 | T_ℓ = d·t_ℓ "forced by the error identity" | ⚠️ | Definitional. Boundedness of t_ℓ is empirical (N ≤ 12). |

## 9. Section 9 — interpretation

| # | Claim | Status | Notes |
|---|---|---|---|
| 9.1 | Overlap vs energy "two faces" | ⚠️ | Interpretive. |
| 9.2 | EE functional Σ H(λᵢ) as entanglement of max-mixed BPS ensemble; "requires Gaussian state" | ⚠️ | H(λᵢ) is *exactly* the single-mode entanglement entropy of decoder eigenvector i. No Gaussianity needed; basis-dependent. LCB parameter pairing reversed (their λ ↔ our a). |
| 9.3 | "Two entropies" S₀ = ln D vs ln d | ⚠️ | Correct observation that exposes the inconsistency in Sec. 3/4 (gravity needs e^{S0} ∝ d). |
| 9.4 | Chaotic vs deterministic: "project out the monotone tower", ⟨r⟩ = 0.601 | ⚠️ | One-flavor model has no monotone tower (that is the two-flavor model). The ⟨r⟩ result itself is CV's (0.599 ± 0.001, N′=12). GUE reference values: surmise 0.60266 vs large-N numerics 0.5996, and should be labelled. |
| 9.5 | Tail is GUE-correlated (⟨r⟩_tail = 0.598 ± 0.010) | ❓ provisional | No code/data in repo. |

## 10. Bibliography of the draft

- `\cite{ChangLin}` arXiv:2412.06902 is Chang–Chen–Sia–Yang, *Fortuity in SYK models*. ❌
- `\cite{FF}` 2608.12160 lacks authors (Chryssanthacopoulos–Vegh). ⚠️
- `\cite{N2matrix}` "Turiaci, Witten et al." — two authors (verified from `papers/turiaci_witten_2023.pdf`). ⚠️
- `\cite{MP}`, `\cite{DistinguishMicro}`: metadata unverified; not in `papers/`. ❓
- SSS, EO: not in `papers/`. ❓ (PSSY, HITZ, BPSmicro = Boruch–Iliesiu–Yan now in `papers/`; see `literature/notes/`.)

## 11. What survives (the "safe core")

1. Exact one-flavor geometry: decoder definition, block BPS equations, endpoint kernels, particle–hole pairing
   (λ ↔ 1−λ at even N′ half filling), the exact norm identities (6.1–6.2).
2. Numerics: GUE adjacent-gap statistics (CV); reproducible, macroscopic, non-Haar redistribution of spectral weight
   toward the walls; atom counts vanish by N = 14.
3. Haar/Jacobi null model: Wachter law, exact Weingarten moments, verified W₁. Clearly labelled as a null model.
4. Kunisky/Chebyshev reparametrization as bookkeeping (not evidence).
