# Kunisky (2023) — Generic MANOVA limit theorems for products of projections

**Citation.** D. Kunisky, *Generic MANOVA limit theorems for products of projections*, arXiv:2301.09543 [math.PR].
File: `papers/kunisky_2023.pdf`.

## Main question
When does the spectrum of ABA, for *structured* (not Haar-rotated) projections A, B, still converge to
Wachter's MANOVA law? What controls the finite-N error?

## Key definitions / equations
- MANOVA(α,β) density with atoms (1 − min(α,β))δ₀ + max(α+β−1,0)δ₁ (eq. 1), edges
  r± = (√(α(1−β)) ± √(β(1−α)))². **Normalization: per ambient dimension N** (includes the zero atom of
  mass 1−α coming from the kernel of A).
- MANOVA(α,β) = Ber(α) ⊠ Ber(β) (eq. 3; Thm. 1.4).
- **Theorem 1.5**: if (1/N)E Tr(A) → α, (1/N)E Tr(B) → β, and (1/N)E Tr[(A−α)(B−β)]^k → 0 for all k, then the
  e.s.d. of ABA converges in moments to MANOVA(α,β) (and A, B are asymptotically free). Thm. 1.7: L² version
  ⇒ convergence in probability.
- Definition 4.2: polynomials q_{k,j,a}(x) (eq. 24).
- Lemma 4.3: **exact algebraic identity for any two projections**: Tr(ABA)^k = β^k Σ_ℓ C(k,ℓ) μ̃_ℓ with a
  recursion for μ̃ (eqs. 25–27) whose sources are N, Tr A − αN, Tr B − βN, Tr[(A−α)(B−β)]^k.
- Lemma 4.4: MANOVA moment recursion (eqs. 46–48); Remark 4.5: explicit moments via Narayana polynomials (see DE15).
- **Lemma 4.6**: Tr(ABA)^k = N·E_MANOVA[X^k] + Δ_k with Δ_k given by the error recursion (eqs. 51–54).
- Thm. 1.12/1.14: largest eigenvalue converges to the edge in the Kesten–McKay case β = 1/2 (moment method
  plus Weingarten bounds).

## Open problems (Sec. 1.4) relevant to us
Problem 5: sufficient conditions for the **smallest non-zero eigenvalue** of non-Jacobi ensembles to converge
to r₋ — exactly the soft-tail question for the decoder.

## Relationship to our project
- The current draft's "exact error formula" (its eqs. errdef–mdev) is built on Lemma 4.6. The citation
  (Lemma 4.6, Definition 4.2) is **accurate**. But:
  1. Lemma 4.3/4.6 are *algebraic identities* valid for any projections. Hence "measured moments = Wachter +
     closed-form error, to machine precision" is a tautology (checked: `scripts/audit/kernel_identity_check.py`
     gives ~1e−14 agreement even for deliberately non-free pairs). It is a change of variables, not a prediction.
  2. Kunisky's MANOVA moments are N-normalized; the current draft's eq. (errdef) writes D·m_k^W with m_k^W
     per-BPS-normalized. The consistent statement is Tr O^k = d·m_k^W + Δ_k (the literal version is off by
     (D−d)m_k^W).
- Theorem 1.5 is stated for **fixed** α, β ∈ (0,1). In SYK, α = a(N) → 0 exponentially, so D-normalized
  mixed traces τ_ℓ = T_ℓ/D → 0 trivially and the theorem says nothing about the per-BPS-state law. The right
  freeness diagnostics are t_ℓ = T_ℓ/d.
- Suggests a clean question: do the SYK mixed traces t_ℓ satisfy t_ℓ → 0 (conditional freeness) or not?
