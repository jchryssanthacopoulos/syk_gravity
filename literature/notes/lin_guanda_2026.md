# Guanda Lin (2026) — Constrained particle on a group: from propagators to correlators

**Citation.** G. Lin, *Constrained particle on a group: from propagators to correlators*, arXiv:2606.30725v1
[hep-th], 29 Jun 2026 (UC Berkeley). File: `papers/lin_guanda_2026.pdf` (84 pp). Journal reference: none at time of
reading.

**Note written 2026-10-05.**
- Read closely: Abstract, §1, §5.2 (four-point functions), §6 (OPE route; Berry curvature).
- Skimmed: §§2–4 and §5.1.
- Not read: Apps. A–D.

## Main question
Can (super-)JT gravity be formulated as a constrained particle on its isometry (super)group, so that supersymmetric
correlators, in particular zero-energy OTOCs, become computable?

## Physical system
N = 0, 2, 4 JT/super-JT. A constrained particle on SL(2,R), SU(1,1|1), PSU(1,1|2). Bulk matter enters as
boundary-anchored Wilson lines.

## Important definitions
- **Physical worldline supercharges.** Bulk supercharges supplemented by compensating large diffeomorphisms (§4.1).
- **Supersymmetric interval propagators** in invariant length variables (§4.2).
- **Three-length composition kernel I₃** (§5.1).
- **Four-point functions.**
  - Time-ordered: A_{4,TOC} = ∫ P₁₂P₂₃P₃₄P₄₁ O(1,2)O(3,4) (5.19).
  - Crossed: A_{4,OTOC} = … O(1,3)O(2,4) (5.20).
  - The diagonals come from hyperbolic geometry (5.21).

## Main assumptions
- Disk level (no wormholes or higher genus).
- Equal-dimension scalar operators for the explicit OTOC.
- Zero energy (BPS) external states for the N = 2 and N = 4 results.

## Main analytical results
1. Super-Schwarzian actions follow from the constrained particle (§3.2). The fermionic constraints explain LMRS's
   auxiliary fermion χ (§3.3).
2. Constraint-preserving supercharges equal the N = 2 / N = 4 super-Poincaré charges (§4.1). Wilson lines are given
   for primaries and descendants (§4.3).
3. Zero-energy 3-pt functions for N = 2 and N = 4 (§5.1).
4. **Four-point functions.**
   - N = 0 OTOC reproduces the JT/Schwarzian OTOC (5.23–5.25; App. D).
   - **N = 2 zero-energy OTOC for equal-Δ scalars is exact:** A₄^{N=2} = Σ_{r=1}^{5} A_r(Δ), with ₂F₁ and ₃F₂ at
     argument −1 (5.26–5.38).
   - Large-Δ saddle at L* = (2−√2)Δ (5.39). The one-loop estimate is about 30 % off at Δ ≈ 1 (Fig. 1).
   - N = 4 BPS: only the large-Δ saddle (5.41). A₄^{N=4}/A₄^{N=2} ∝ Δ⁶ (5.42).
5. **§6, OPE route.** The TOC as a sum over channel representations (6.1). The OTOC needs a further exchange matrix
   Λ_Ξ (6.8). The super-6j / exchange data are not yet known.
6. **§6, Berry curvature.**
   - F_IJ = ⟨a|[∂_IΠ_BPS, ∂_JΠ_BPS]|b⟩ (6.9), or its spectral form (6.10).
   - The second moment V_{IJ;IJ} = (1/D_BPS) Tr F_IJ F_IJ† (6.15) equals the difference of two orderings, IJJI − IJIJ
     (6.16). In gravity this is TOC − OTOC with BPS and non-BPS propagators and energy denominators (6.17–6.19).
   - Proposed for D1–D5–P with the marginal operator as a superdescendant (6.20).

## Important equations for us
(5.19)–(5.22), (5.26)–(5.38), (6.9)–(6.10), (6.15)–(6.19).

## Numerical methods
Evaluation of the exact hypergeometric sums against the large-Δ saddle (Fig. 1). No SYK comparison.

## Limitations
- Disk-level super-JT, i.e. the conformal (Schwarzian) regime.
- At our sizes (N′ ≤ 16) the super-Schwarzian overshoots the measured two-point function of our probe by 15–60 %,
  whereas the double-scaled chord theory (BLY) matches it to 1–2 % (see BLY note). Quantitative use at our N′ is
  therefore not expected.
- Which R-charge sector and operator type (neutral vs charged scalar) the N = 2 OTOC formula refers to must be
  checked in §4–5 before any comparison.
- The normalization relative to the TOC is not stated as a ratio.

## Relationship to our project
1. **This is the first explicit zero-energy N = 2 crossed four-point function in our collection.**
   - LMRS give only the large-Δ ratio OTOC/TOC ~ e^{−3.52Δ} (their eq. 159).
   - BLY set up the one-particle wormhole but do not evaluate it (their §5).
   - Lin's (5.33)–(5.38) at Δ = 1/3 is therefore the super-JT prediction for the crossed/uncrossed ratio of our
     light probe, once normalized by the TOC (5.19). It is useful as the triple-scaled limit of whatever chord-level
     answer we obtain for the decoder remainder R.
2. **Berry curvature: the analogy is structural, not literal** [derived in project, `docs/derivations.md` D5].
   - Our OU fourth moment expanded about x = 1 is
     T₄ − 2T₂ + a = τ((PP_xP − P)²) ≈ (4ε²/D) E[Tr g_tot² + 2Σ_IJ Tr g_IJ²], where ε = 1 − x and g is the
     operator-valued quantum metric.
   - It contains **only the quantum metric**. The IJIJ and IJJI orderings enter with equal weight, and the
     Berry-curvature pieces cancel.
   - So isotropic OU noise does *not* measure Lin's TOC − OTOC combination. A curvature measurement needs two distinct
     deformation directions applied in a loop.
3. **The decoder is a conjugation deformation.**
   - P′ = e^{iπn_i} P e^{−iπn_i} is the θ = π endpoint of the Chen–Colin-Ellerin–Mamroud–Papadodimas conjugation
     flow with Λ = n_i.
   - Along that flow the quantum metric is Λ̂(1 − Λ̂) with Λ̂ = P n_i P = P − O_T.
   - Lin's §6 protocol therefore applies to *pairs* of sites (Λ_μ = n_i, Λ_ν = n_j). A single-site decoder is a
     one-direction (metric) quantity.
