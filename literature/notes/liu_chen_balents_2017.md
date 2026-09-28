# Liu, Chen & Balents (2017) — Quantum entanglement of the Sachdev–Ye–Kitaev models

**Citation.** C. Liu, X. Chen, L. Balents, *Quantum Entanglement of the Sachdev-Ye-Kitaev Models*,
arXiv:1709.06259 (v3, 2018); Phys. Rev. B 97, 245126 (2018) as cited in the project draft (journal ref not
verified from the PDF). File: `papers/liu_chen_balents_2017.pdf`. (Read: Secs. I–III.)

## Main results
- For free SYK₂, the subsystem correlation matrix C_A = V†V (V a k×m block of a Haar unitary) is β=2 Jacobi,
  and its eigenvalue density is the Wachter law
  f(x,κ,λ) = (1/(2πλ)) √((λ₊−x)(x−λ₋))/(x(1−x)) + (1−κ/λ)Θ(λ−κ)δ(x), λ = m/N, κ = k/N,
  λ± = (√(κ(1−λ)) ± √(λ(1−κ)))² (eq. 2).
- Entanglement entropy S_A = m ∫ [−x ln x − (1−x) ln(1−x)] f dx (eq. 3), valid because the state is Gaussian.
- SYK₄ ground-state EE is below the Page value when subsystems are comparable; crossover SYK₂ → SYK₄ with energy.

## Relationship to our project
- Cited by the current draft as an "open bridge". Two corrections:
  1. **Parameter mapping:** LCB normalize by λ = m/N, the dimension of the *compressed* space. So LCB's λ ↔ our a
     (BPS fraction) and LCB's κ ↔ our b (slot fraction). The draft states the reverse pairing "(κ,λ) in place
     of (a,b)".
  2. The Gaussian-state caveat is misplaced for the decoder. For a decoder eigenvector Ψ = a + b∧e, the reduced
     state of the single added mode is diag(λ, 1−λ) (cross terms vanish by fermion-parity superselection). So
     H(λ_i) is *exactly* the entanglement entropy of eigenvector i across the (old modes | new mode) cut.
     Σ_i H(λ_i) is therefore an exact average single-mode entanglement of the decoder eigenbasis, not a
     free-fermion proxy. (It is basis-dependent: generic BPS states have entropy H(⟨Π⟩).)
