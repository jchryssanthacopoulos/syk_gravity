# Liu & Shen (2026) — Limit joint distributions of SYK models with partial interactions, mixed q-Gaussians, ε-freeness

**Citation.** W. Liu, H. Shen, *Limit joint distributions of SYK Models with partial interactions, Mixed q-Gaussian
Models and Asymptotic ε-freeness*, arXiv:2602.00789v2 [math.OA], 5 Apr 2026 (v1 31 Jan 2026). File:
`papers/liu_shen_2026.pdf` (29 pp). Read on 2026-09-30: §§1–4 in full, the §5 proofs skimmed.

## Main question
What is the large-n joint distribution of several SYK Hamiltonians with *different* interaction lengths, supported on
*partially overlapping* sets of Majorana fermions? Can such models realize every mixed q-Gaussian system, and in
particular ε-freeness (graph products)?

## Physical system
Majorana SYK (Feng–Tian–Wei normalization):
- H_{k,n} = (√−1)^{⌊r/2⌋} C(n,r)^{−1/2} Σ_{I⊂A_k} J_{k;I} Ψ_I;
- the couplings are independent across k and I, with mean 0, variance 1 and bounded moments;
- |A_k| = n;
- interaction lengths r_{k,n} with r_{k,n}/n → 0 and fixed parity.

## Important definitions
- **Mixed q-Gaussian (Q-)system:** l_i* l_j − q_ij l_j l_i* = δ_ij, s_i = l_i + l_i* (Bożejko–Speicher; Speicher).
  Moments are τ(s_{ε(1)}…s_{ε(d)}) = Σ_{π∈P₂, π≤ker ε} Π q_ij^{cr(π;i,j)} (Prop. 2.3), where cr counts crossings
  between an i-block and a j-block.
- **ε-freeness (Mlotkowski; Def. 3.1).** A mixture of classical (ε_ij = 1: commuting) and free (ε_ij = 0)
  independence, realized by graph products. It coincides with a Q-system with q_ij ∈ {0, 1} (Cor. 3.4, Prop. 3.5).

## Main results
- **Thm. 1.1** (independent SYK on the same fermions). If r_i r_j/n → λ_ij, the family converges to a Q-system with
  **q_ij = (−1)^{r_i r_j} e^{−2λ_ij}**. This is Pluma–Speicher for equal lengths, extended to different lengths.
- **Thm. 1.2** (partial overlaps). With supports A_i, A_j and r_i r_j |A_i ∩ A_j|/n² → λ_ij, the same formula holds.
  Varying the overlap realizes any q_ij ∈ [0, 1]. Negative q_ij needs odd lengths, with the sign constraint noted in
  §5.4.
- **§3 end.** An explicit construction of **asymptotic ε-freeness** (q_ij ∈ {0,1}): take disjoint supports for
  q_ij = 1, and overlaps with r_i r_j |A_i∩A_j|/n² → ∞ for q_ij = 0. This answers an open question of
  Morampudi–Laumann.
- **Key lemmas.**
  - Lemma 4.4: if any λ_ij = ∞, the sign average E(−1)^{Σ|R_i ∩ R_j|} → 0.
  - Lemma 4.5: the factorial moments of the overlaps are asymptotically Poisson, giving
    E(−1)^{Σ|R_i∩R_j|} → Π e^{−2λ_ij}.
- **Example 4.6.** If r/n ↛ 0 (e.g. r = n/2), the Poisson limit and Thm. 1.2 **fail**: E(−1)^{|R₁∩R₂|} = 1/2 vs the
  predicted e^{−1/2}.
- **Open questions (§5.4):** (1) the q_ij relations when some r/n ↛ 0; (2) realizing all mixed Q-systems.

## Numerical methods
None; the proofs are combinatorial (moment method).

## Limitations
- Moments of the *Hamiltonians* only: no projectors, no ground states, no finite-n corrections.
- Couplings are independent across the family; there are no correlated couplings.
- Only Majorana, non-supersymmetric models.

## Relationship to our project
1. **Rigorous version of the "inter-family crossing weight" picture.**
   - Two independent SYK Hamiltonians on the *same* fermions with the same p have **q₁₂ = q₁₁ = e^{−λ}** (Thm. 1.1):
     they cross each other exactly as they cross themselves.
   - Freeness (q₁₂ = 0) needs λ₁₂ → ∞.
   - So at finite λ two independent SYK models are *not* free, the Hamiltonian-level counterpart of our finding that
     two independent SYK BPS spaces are not mutually free.
   - It also explains the **trend** [interp]. At fixed p = 3, λ = 18/N′ falls as N′ grows, so q₁₂ → 1 (towards
     classical/commuting). The two BPS spaces should become *less* free with N′, as observed: E τ(P₁P₂P₁P₂) is
     +0.6, +2.0, +4.3 % above free at N′ = 10, 12, 14.
2. **Lemmas 4.4–4.6 are the rigorous content of our finite-N caveat.**
   - The Z₂ weight x per joining chord is an average (−1)^{|I∩S|}. Factorization over several chords (x^{#chords})
     is exactly their Poisson limit (Lemma 4.5).
   - Our parity probe at s = N′/2 is their Example 4.6 regime (|S|/N′ ↛ 0). There the exact per-chord weight is 0,
     while the double-scaled Poisson value is e^{−2ps/N′} = e^{−3} ≈ 0.05 at N′ = 14 (`parity_family_table.md`:
     x_exact = 0, x_ds = 0.050). The failure of the Poisson limit is visible directly in our data.
3. **What they do not cover, and we need.**
   - *Correlated* couplings: Q′ = UQU has the *same* C_I with signs. In Q-system language, Q and Q′ are q-Gaussians
     with nonzero covariance x, not independent generators with a crossing parameter.
   - Supersymmetric structure, and BPS projectors/compressions rather than Hamiltonian moments.
   - A mixed-q-Gaussian extension with correlated generators would be the natural mathematical home for the
     decoder's four-arc diagrams (research question on R).
4. Useful for citations: Pluma–Speicher (independent SYK → q-Gaussian systems) and Feng–Tian–Wei (SYK spectrum is
   q-Gaussian, q = (−1)^{r} e^{−2λ}) are their refs. [27] and [14].
