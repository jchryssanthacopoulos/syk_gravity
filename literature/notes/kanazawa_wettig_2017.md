# Kanazawa & Wettig (2017) — Complete random matrix classification of SYK models with N = 0, 1 and 2 supersymmetry

**Citation.** T. Kanazawa, T. Wettig, *Complete random matrix classification of SYK models with N = 0, 1 and 2
supersymmetry*, arXiv:1706.03044 (v2, 30 Aug 2017) [hep-th], RIKEN-QHP-311. The PDF says "prepared for submission to
JHEP"; the journal reference has not been checked. File: `papers/kanazawa_wettig_2017.pdf` (42 pp).

**Note written 2026-10-26.**
- Read in full: §§6, 7 and Apps. A–B.
- Skimmed: §§1–5 (N = 0, 1 SYK and the hybrid model of §5).
- CV cite this paper as their ref. [25], for the name "exceptional" of the window-boundary BPS states.

## Main question
What are the Altland–Zirnbauer symmetry classes of SYK with N = 0, 1, 2 SUSY? In particular, for N = 2, what are the
level statistics in each fermion-number eigenspace, and what explains the number of zero modes?

## Physical system
- N=2 SYK: N_c complex fermions, q̂ = 3 (and general odd q̂ in §6.5), H = {Q, Q̄}.
- Q is a q̂-body supercharge with independent complex Gaussian couplings (eq. 6.1).
- Its orientation matches ours: Q raises the fermion number f by 3, as our Q = ω∧ does (checked below).
- §§3–5 treat N = 0 and N = 1 SYK and a new non-SUSY "hybrid" model.

## Important definitions
- **Mod-6 blocks.** V_f is the fermion-number eigenspace, of dimension C(N_c, f).
  - M = Q + Q̄ preserves f mod 3 and anticommutes with (−1)^F. This gives three chiral blocks V^f ⊕ V^{f+3}.
  - The block dimensions D_f (f mod 6) are given in eq. 6.4 and App. A.
- **Three-way split of each V_f** (eqs. 6.8–6.10): V_f = V_f⁺ ⊕ V_f⁻ ⊕ V_f^z, with dimensions N_f⁺, N_f⁻, N_f^z.
  - V_f⁺: Qψ ≠ 0, Q̄ψ = 0. **In our language: the coexact forms of degree f (im Q†).**
  - V_f⁻: Qψ = 0, Q̄ψ ≠ 0. **In our language: the exact forms (im Q).**
  - V_f^z: zero modes. **In our language: the harmonic (BPS) forms.**
  - Hence N_f⁺ = rank(Q: Λ^f → Λ^{f+3}), N_f⁻ = rank(Q: Λ^{f−3} → Λ^f), and N_f^z = dim H^f(Q).
- **Relations.**
  - Particle–hole P gives N_f⁺ = N⁻_{N_c−f} and N_f^z = N^z_{N_c−f} (eq. 6.11).
  - Q maps V_f⁺ onto V⁻_{f+3} bijectively, so N_f⁺ = N⁻_{f+3} (eq. 6.12). This is the coexact–exact pairing of
    supermultiplets.

## Main assumptions
- **Premise (6.15).** All zero modes lie in |f − N_c/2| ≤ 3/2, with equality only for "exceptional" zero modes when
  N_c ≡ 1 (mod 4). This is CV's fortuity window (CV eq. 12; CCSY).
  - The support is numerical, from Schwarzian arguments, and from a heuristic "pinching" argument (§6.4).
  - The pinching argument: along an exact sequence ⋯ → V_f → V_{f+3} → ⋯, a generic nilpotent map has maximal rank, so
    cohomology survives only at the largest-dimension term.
  - This is the generic-position logic of CV App. A and of our counting laws (D8 §6).
- **Generic couplings.** All tables hold "almost surely" (footnote 16).

## Main analytical results
1. **Zero-mode bound** (eq. 6.5). N_z ≥ Σ_{f=0,1,2} |D_f − D_{f+3}|, which is 4·3^{N_c/2−1} for even N_c and
   2·3^{(N_c−1)/2} for odd N_c. This equals the Witten-index bound. It is saturated for N_c ≡ 0, 2, 3 (mod 4).
2. **Exceptional zero modes explained** (§6.2).
   - For N_c ≡ 1 (mod 4), the block V^f ⊕ V^{f+3} that contains f = (N_c ∓ 3)/2 is in class DIII (β = 4). Kramers
     degeneracy forces every level of M to be doubly degenerate.
   - In that block D_f = D_{f+3} = (2^{N_c−1} − 1)/3, which is odd. Kramers pairs therefore cannot fill it, so there are
     at least 2 zero modes, at f = (N_c ± 3)/2.
   - Footnote 12: N_c = 5, 13, 17 have 2 exceptional zero modes, but N_c = 9 has 6. "The origin of the 4 additional
     zero modes is unclear."
3. **Exact sequences** (Table 8). Given (6.15), each residue class f mod 3 gives a sequence V_{f₀} → V_{f₀+3} → ⋯
   that is exact except at its window term.
   - Closed forms follow: N_f⁺ = Σ_n (−1)^n C(N_c, f − 3n) below the window (eqs. 6.16–6.17), which is our D̃_p (CV
     eq. 114).
   - Near the window: N_f⁺ = C(N_c, f) − N⁺_{f−3} − N^z_f (eqs. 6.20–6.21).
4. **BPS counts per sector** (eq. 6.19).
   - N_c ≡ 0, 2 (mod 4): 2·3^{N_c/2−1} at f = N_c/2 and 3^{N_c/2−1} at N_c/2 ± 1.
   - N_c ≡ 1 (mod 4): 3^{(N_c−1)/2} at (N_c ± 1)/2, plus 1 or 3 at (N_c ± 3)/2.
   - N_c ≡ 3 (mod 4): 3^{(N_c−1)/2} at (N_c ± 1)/2.
5. **RMT classes** (eq. 6.13). H must be analysed on V_f⁺ and V_f⁻ separately, because QQ̄ and Q̄Q commute with H.
   - For even N_c every V_f^± is GUE.
   - For odd N_c the sectors V⁺_{(N_c−3)/2} and V⁻_{(N_c+3)/2} are special: GOE for N_c ≡ 3 and GSE for N_c ≡ 1
     (mod 4). The antiunitary operator PQ maps them to themselves. Every other sector is GUE.
   - These two sectors form **one supermultiplet**, the one that straddles the window from f = (N_c−3)/2 to (N_c+3)/2.
   - Positive levels are 4-fold degenerate (quadruplets ψ, Qψ, Pψ, PQψ; eq. 6.14). At even N_c, H on V_{N_c/2} is two
     identical GUE copies (footnote 13).
6. **General odd q̂** (§6.5, eq. 6.24). The special sectors move to f = (N_c ∓ q̂)/2. Zero modes lie in
   |f − N_c/2| ≤ q̂/2, with equality only at N_c = 7, 11 (2 modes each) for q̂ = 5.

## Important equations
6.4–6.5, 6.8–6.13, 6.15–6.21, 6.24, App. A (A.1–A.6), App. B tables.

## Numerical methods
Exact diagonalization, using ratio-of-spacings distributions (Atas et al.). 10²–10³ realizations at N_c = 15–17
(Figs. 10, 12). App. B tabulates N_f^±, N_f^z for N_c = 3–17 at q̂ = 3.

## Relevant figures / tables
- Fig. 11: the exact-sequence structure at N_c = 15, with dimensions.
- Table 8: the exact sequences by N_c mod 6.
- App. B: per-degree Hodge dimensions.

## Limitations
- Premise (6.15) is not proven.
- The N_c = 9 excess of exceptional modes is unexplained.
- The classification uses short-range level statistics only.
- Nothing is said about eigenvectors, operators projected to the BPS space, or how BPS spaces at different N_c are
  related.

## Notation warning: CV App. A uses the opposite ± labels
- CV App. A (eqs. 110–111) defines V⁻_p = im Q† ⊕ B_e,− (coexact, plus the exceptional states at ℓ = −q) and
  V⁺_p = im Q ⊕ B_e,+ (exact, plus those at ℓ = +q). KW's V⁺_f is the coexact space and V⁻_f the exact one.
- So CV's V^∓ correspond to KW's V^± (KW's tables count the exceptional states separately, in N^z).
- CV's Fig. 15 statement that its walk counts T_p (exact) and B_p (coexact) "match those in Appendix B of [25]" holds
  with this swap. Checked against the N_c = 5 row: KW N⁺ = (1, 4, 1, 0, 0, 0) = rank(Q out of degree f) = coexact.
- Miyahara–Shibuya 2026 (their eqs. 7–9 and ref. [21]) use KW's three-way decomposition and table (6.13) to choose the
  GUE sector V⁺_6 at N = 12.

## Relationship to our project
1. **Hodge dimensions for every degree.**
   - App. B gives exactly the exact, coexact and harmonic dimensions that the uplift analysis needs at every degree, not
     only at half filling.
   - [verified 2026-10-26] Our code (`chord_chaos.WedgeStructure`, one random realization) reproduces the App. B rows
     for N_c = 5, 9, 12, 13 exactly, under the dictionary above: V⁺ = coexact, V⁻ = exact, V^z = BPS. This includes the
     exceptional states (1 at N_c = 5, 13; 3 at N_c = 9).
   - Our `uplift.exact_forms_dim` / `coexact_forms_dim` recursions are KW eq. (6.16) and its particle–hole mirror. They
     are valid outside the window.
2. **Exceptional states.**
   - CV's name and their explanation (CV eqs. 14–15: an odd-dimensional antisymmetric map M) follow KW §6.2 (the DIII
     Kramers argument).
   - The off-half-filling uplift analysis (D10 §3) gives a third, inductive version of the same fact.
     - The exceptional states at N′ are the radical of the connecting form B(β, α) = ∫β∧η∧α on the ordinary BPS
       space H(N′−1, 1−q).
     - B has symmetry (−1)^k with k = (N′−q)/2, so a radical is forced when k is odd and dim H is odd.
     - For q = 3 this means N′ ≡ 1 (mod 4) on the odd dimension 3^{N/2−1}.
   - [verified 2026-10-26]
     - Radicals 1, 3, 1 at N = 4, 8, 12, and none at N = 6, 10, 14 (q = 3).
     - q = 5: radical 1 at N = 6, 10. At N = 14 the form is antisymmetric but the dimension is even (1000), so there is
       no radical, in agreement with KW's "none at N_c = 15".
   - **Beyond KW's tested range (q̂ = 5).**
     - The parity criterion forces exceptional states again at N_c = 19, 23, 31, 35, 43, …, and none at 15, 27, 39.
     - For N_c = 19 the statement is rigorous for generic couplings: CV's M on Λ⁷(19) is antisymmetric on the even
       dimension 50388, and Q is injective on Λ² (rank 171, checked for two realizations). So dim H⁷ ≡ 171 ≡ 1 (mod 2).
     - KW's "saturated only for N_c = 7 and 11" holds only up to their N_c = 17.
   - KW's N_c = 9 excess is the rank deficiency of B on H(8, −2): radical 3 instead of 1.
     - The radical moves with η, but its dimension is 3 for every η tried.
     - It remains unexplained.
3. **The special supermultiplet is the window-edge multiplet.**
   - KW's GOE/GSE pair V⁺_{(N_c−3)/2} ↔ V⁻_{(N_c+3)/2} is the coexact family at the lower window edge and the exact
     family at the upper edge.
   - These are the edges s = ±1/2 of the chord window (BBNR; D9), where the chord continuum is gapless (HITZ App. A: the
     Q = ½ multiplet).
   - [derived, D10 §4] The first-order dressing of the uplift into the off-centre sectors |ℓ′| = q − 1 passes through
     exactly this multiplet. K = Q₀⁺E takes the ℓ = ∓1 slot to the exact (coexact) forms at ℓ = ±q, and back to the
     other edge.
   - [num] Its Hodge spectrum has a hard edge at zero: E_min/E_mean = 4e−6 at N = 11 and 9e−7 at N = 13. The dressing
     cost in those sectors therefore has a heavy tail. Note also that Q₀ maps between spaces of equal dimension only
     for this multiplet: C(N, (N−q)/2) = C(N, (N+q)/2).
4. **Decoder statistics.**
   - KW's "every V_f^± is GUE at even N_c" is consistent with CV's ⟨r⟩ = GUE for the decoder at N′ = 12, P = 6.
   - It also bears on the Turiaci–Witten note's caution about β = 1 in particle–hole-symmetric sectors. That exception
     belongs to odd N_c, i.e. to the edge multiplet, not to half filling at even N_c.
   - At even N_c, H on V_{N_c/2} is two identical GUE copies (P pairs V⁺ with V⁻). This is the same particle–hole
     structure that gives the exact λ ↔ 1−λ pairing of the half-filling decoder at even N′ (D6 §5).
