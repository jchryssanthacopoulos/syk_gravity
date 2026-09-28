# Chen, Lin & Shenker (2024) — BPS Chaos

**Citation.** Y. Chen, H. W. Lin, S. H. Shenker, *BPS Chaos*, SciPost Phys. 18, 072 (2025),
arXiv:2407.19387 (v4, May 2026). File: `papers/chen_lin_shenker_2024.pdf`.
(Reading depth: Secs. 1, 5, 7 read in full; Secs. 3, 4, 6 skimmed.)

## Main question
Can the LMRS diagnostic separate BPS sectors described by horizonless geometries (fuzzballs, LLM, LM)
from those described by macroscopic BPS black holes?

## Physical system
1/2- and 1/4-BPS sectors of N=4 SYM at weak coupling (one-loop dilatation operator D₂); 2-charge
(D1-D5) sector of Sym^N CFTs; qualitative discussion of 1/16-BPS.

## Key definitions
- **LMRS operator** Ô = P_BPS O P_BPS (eq. 1.2); **LMRS criterion**: chaotic BPS subspace ⇔ Ô shows RMT
  behaviour. **Strong** chaos = Thouless time O(1); **weak** = t_Th growing as a power of N.
- ETH form ⟨E_i|O|E_j⟩ = f(E)δ_ij + e^{−S/2}R_ij (eq. 1.3) applied to a degenerate subspace.

## Main results
1. In an N=1 SUSY random-matrix ensemble, BPS states are the last ν columns of a Haar unitary, hence a
   **random hyperplane** when ν ≪ L (eqs. 1.4–1.5); footnote 8: when ν is comparable to L, correlations
   among the vectors matter.
2. From super-Schwarzian: width of the Ô_Δ spectrum ~ 1/S0^Δ because ⟨ÔÔ⟩ ~ e^{−Δℓ} with typical
   wormhole length ℓ ~ 2 log S0 (discussion around eq. 1.8); for Δ ≫ 1 a semicircle, otherwise a more
   general RMT potential (footnote 11).
3. Conjecture (Sec. 1.3): horizonful BPS sectors show strong LMRS chaos; horizonless ones do not.
4. 1/2-BPS (free fermions) and 1/4-BPS numerics: projected operators are banded ⇒ weak chaos with t_Th ~ N^#.
5. Sec. 5: fortuity "invades" the BPS subspace with chaotic non-BPS states as N decreases (qualitative).
6. Sec. 6: projected twist operator in D1-D5 (eqs. 6.27–6.29) — v4 note: a later complete treatment removed the
   weak chaos found here.

## Limitations
No first-principles strong-chaos computation in any boundary theory; 1/16-BPS argument is indirect.

## Relationship to our project
- Supplies the *definition* of "BPS chaos" behind CV's ⟨r⟩ = GUE result for O_T.
- The ν ≪ L caveat is directly relevant: our BPS fraction a = d/D ≈ 0.3–0.6 at accessible N, so a
  "random hyperplane" (Haar/Wachter) model is *not* justified by this argument at those sizes.
- The 1/S0^Δ width statement supports the hypothesis that the physical decoder variance decays as a power of
  N, not exponentially like the Wachter variance (see LMRS note).
- Thouless-time analysis of O_T has not been done in our project.
