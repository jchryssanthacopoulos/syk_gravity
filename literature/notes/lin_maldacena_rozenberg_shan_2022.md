# Lin, Maldacena, Rozenberg & Shan (2022) — Looking at supersymmetric black holes for a very long time

**Citation.** H. W. Lin, J. Maldacena, L. Rozenberg, J. Shan, *Looking at supersymmetric black holes for a
very long time*, SciPost Phys. 14, 128 (2023), arXiv:2207.00408 (v5, Aug 2024). Companion: *Holography for
people with no time*, arXiv:2207.00407 (not in `papers/`). File: `papers/lin_maldacena_rozenberg_shan_2022.pdf`.

**Role in this project.** The origin of "LMRS operators" Ô = POP. Our decoder O_T = P_B Π_T P_B is of this
type (Π_T = 1 − n_{N′}). LMRS supply the *actual* N=2 super-JT computation of BPS-projected correlators — the
calculation any gravitational account of the decoder must reproduce.

## Main question
What do simple operators look like after projection onto the zero-energy (BPS) subspace of an extremal
N=2 black hole, and how can their correlators be computed in N=2 super-JT / super-Schwarzian?

## Physical system
N=2 super-Schwarzian (boundary mode of N=2 super-JT) with bulk matter; tested against N=2 SYK (q̂=3) at N=16
complex fermions.

## Key definitions
- Ô ≡ P O P, P ≡ lim_{β→∞} e^{−βH}; observables Tr[Ô₁⋯Ôₙ] (eq. 1).
- Two-sided super-Liouville QM for the wormhole length ℓ and R-phase a (eqs. 4–10); zero-energy bound
  states |Z_j⟩, |j| < 1/2 (eq. 22), norm 1/cos πj (eq. 26).
- Disk partition function Z(β) with BPS part e^{S0} Σ_{|j|<1/2} cos πj (eq. 38).
- N_BPS = Tr P = e^{S0} L̂, L̂ = Σ_{|j|<1/2} cos πj (eq. 82).

## Main analytical results
1. Zero-energy two-point functions (the long-time constants):
   charged BPS operator ⟨Ô†Ô⟩_j (eq. 65); **neutral** operator
   ⟨ÔÔ⟩_j = (e^{S0}/π)(1/2)(cos πj)² ΔΓ(Δ)² Γ(Δ+1/2±j)/Γ(2Δ) (eq. 66).
2. SYK normalization: projected correlators carry (2α_S N)^{−2Δ} (eqs. 83–85; α_S = 0.00842 for q̂=3,
   eq. 80). E.g. the neutral bilinear ψ^kψ̄^i (Δ = 2/6) in eq. (85).
3. Large-Δ: OTOC/TOC ~ exp(−3.52…Δ) (eq. 159) — **crossed diagrams are suppressed by Δ**, not by e^{S0}.
   Then Ô is Gaussian L×L with a semicircle (Sec. 5.3, eqs. 163–166; Rényis give Catalan numbers, eq. 166);
   entanglement of a one-particle wormhole S = log L − 1/2 (eq. 165).
4. Δ → 0: every Wick contraction (crossed or not) has equal weight: Z_n/Z_1^n = L^{1−n}(2n−1)!! (eq. 169),
   S − S0 = −0.7296….
5. Sec. 5.4: for Δ ≫ 1 the connected SFF of Ô has a linear ramp; expected by universality for Δ > 0.
6. BPS operators placed a long time apart are "on top of each other" (eqs. 86, 161).

## Numerical methods / results
ED of N=2 SYK, N=16: Table 1 compares predicted vs measured zero-energy two-point functions per R-charge
(e.g. neutral ψ̄_iψ_j: 0.0874 predicted vs 0.079 ± 0.001). Fig. 12: spectrum of P ψ₁ψ̄₂ P at N=16 is dense
but "not quite a semi-circle"; extra zero modes at small N (App. G; absent for N ≥ 20).

## Limitations
Cylinder results (Sec. 2.7) neglect particles wrapping the cylinder; reliable only for Δ ≫ 1. Finite-Δ
spectra of Ô are not computed. SYK comparison is at a single N.

## Relationship to our project
- **Contradicts two load-bearing claims of the current draft.** (i) "Crossing chords = handles ∝ 1/D²":
  in LMRS, crossings on the *disk* are weighted by Δ-dependent OTOC factors (O(1) in e^{S0}); only genus is
  suppressed by e^{−2S0}. (ii) "BPS projection trivializes the Schwarzian so only e^{S0χ} survives": the
  zero-energy correlators (65)–(66) are nontrivial functions of Δ and j.
- **Entropy identification:** e^{S0} counts BPS states (eq. 82), i.e. corresponds to our d, not to the
  ambient Fock-sector dimension D used as e^{S0} in the current draft.
- n_{N′} − ⟨n⟩ is a neutral fermion bilinear with Δ = 2 × 1/(2q̂) = 1/3 for q̂=3. Eq. (85)-type formulas
  therefore predict the **variance of the decoder spectrum** m₂ − m₁² and its N-dependence ~ N^{−2Δ}
  (power law), in contrast with the Haar/Wachter variance ab(1−b) ∝ a ~ (√3/2)^N (exponential). This is the
  top open calculation (see `docs/research_questions.md`, Q1). Normalization details (which j-sectors,
  i = k vs i ≠ k, odd N′) still need to be worked out.
- Our Π_T is a light operator (Δ = 1/3), so LMRS's Δ→∞ semicircle regime does not apply directly.
