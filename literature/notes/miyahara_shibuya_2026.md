# Miyahara & Shibuya (2026) — Chaos–integrability transition in the BPS subspace of the N=2 SYK model

**Citation.** L. Miyahara, S. Shibuya, *Chaos-Integrability Transition in the BPS Subspace of the N=2 SYK
Model*, arXiv:2605.20913 (v2, Aug 2026). File: `papers/miyahara_shibuya_2026.pdf`.

## Main question
Can a chaos→integrability crossover be diagnosed *purely inside the BPS subspace*?

## Physical system
Q_g = Q_com + 2^{−gN/2} Q_SYK (eq. 15), where Q_SYK = Σ C_ijk ψ_iψ_jψ_k (eq. 2) and
Q_com = Σ_a C_a ψ_{3a−2}ψ_{3a−1}ψ_{3a} (eq. 13) gives a commuting, integrable H_com (eq. 14).
N = 12, sector f = 6.

## Key points
- Hilbert space decomposition H_f = H_f⁺ ⊕ H_f⁻ ⊕ H_f^BPS (eq. 7); non-BPS sectors are GUE for Q_SYK (Figs. 1–2,
  ⟨r⟩ ≈ 0.598).
- Probe O = Σ_i A_i ψ_i†ψ_i with random real A_i (a random combination of occupation numbers), projected to BPS.
- ⟨r⟩, P(s) and Gaussian-filtered SFF of O_BPS cross over from GUE to Poisson around g ≈ 1 (Figs. 5–7).
- Both Q_SYK and Q_com have only **fortuitous** BPS states, yet Q_com is integrable.

## Relationship to our project
- Uses essentially our kind of operator (occupation numbers projected to BPS).
- Shows that "fortuitous" does **not** imply "chaotic". This directly tests CV's conjecture that metric
  fragility signals BPS chaos: the decoder should be computed along the Q_g interpolation (open question).
- The random weights A_i were chosen to avoid degeneracy; our single-mode Π_T is non-degenerate on the BPS
  space only generically — worth checking at the integrable end.
