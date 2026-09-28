# Collins, Matsumoto & Novak (2021) — The Weingarten calculus

**Citation.** B. Collins, S. Matsumoto, J. Novak, *The Weingarten Calculus*, arXiv:2109.14890 [math-ph]
(Notices-style survey). File: `papers/collins_matsumoto_novak_2021.pdf`. (Read in full; Secs. 5–6 lightly.)

## Main content
- Fundamental theorem of Weingarten calculus (Thm. 2.1): Haar integrals of matrix-element monomials are
  expressed through the Gram matrix of invariants and its (pseudo)inverse, the Weingarten matrix.
- Unitary group (Sec. 4): I_ij = Σ_{ρ,σ∈S(d)} δ_{i,i′ρ} δ_{j,j′σ} W_ρσ (Thm. 4.4). Character formula via Jucys–Murphy
  elements (Thm. 4.3).
- **1/N expansion (Thm. 4.5):** for 1 ≤ d ≤ N,
  W_ρσ = (−1)^{|ρ⁻¹σ|} N^{−d−|ρ⁻¹σ|} Σ_{k≥0} W̃_k(ρ,σ) N^{−2k},
  where W̃_k counts **weakly monotone walks** on the Cayley graph of S(d) from ρ to σ of length |ρ⁻¹σ|+2k.
  Monotone walks are the "Feynman diagrams" of Haar integrals. Leading term W̃₀ is a product of Catalan numbers
  (Thm. 4.6), i.e. the Möbius function of non-crossing partitions.
- Orthogonal/symplectic analogues with pairings and hyperoctahedral groups (Sec. 5); applications including
  asymptotic freeness (Sec. 6).

## Relationship to our project
- The correct reference for the Haar null model's genus expansion (draft §2, verified in
  `scripts/audit/weingarten_W1_check.py`, which inverts the Gram matrix exactly as in Thm. 2.1).
- Thm. 4.5 makes precise the audit's objection. The 1/D² corrections of the Haar model are **signed counts of
  monotone walks** (alternating signs from (−1)^{|ρ⁻¹σ|}). They are not a positive sum over surfaces weighted by
  D^χ, as the draft's "fatgraph duality with the super-JT action" asserts. The negative weight −ab² in the draft's
  "four surfaces of m₂" is exactly such a Möbius/Weingarten sign.
- Thm. 4.6 / Catalan leading term is the "Fuss–Catalan fingerprint" the draft mentions; it is a property of the
  Haar model, not evidence for gravity.
