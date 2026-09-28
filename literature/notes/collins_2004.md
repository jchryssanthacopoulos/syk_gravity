# Collins (2004) — Product of random projections, Jacobi ensembles and universality problems

**Citation.** B. Collins, *Product of random projections, Jacobi ensembles and universality problems arising
from free probability*, arXiv:math/0406560 [math.PR]. File: `papers/collins_2004.pdf`.
(Journal version not verified from the PDF.)

## Main results
- **Theorem 2.2 (exact, finite n)**: for a fixed projection π of rank q_n and a Haar-random projection π̃ of rank
  q̃_n, the compression ππ̃π (on ran π) has the law of J = (X+X′)^{−1/2} X (X+X′)^{−1/2} with independent
  complex Wisharts; for q̃_n ≥ q_n, q_n + q̃_n ≤ n it is a Jacobi unitary ensemble of parameters
  (q_n, n − q_n − q̃_n, q̃_n − q_n).
- Sec. 3.2: asymptotic freeness ⇒ limit µ₁ ⊠ µ₂ with the Wachter density
  √((r₊−x)(x−r₋))/(2πx(1−x)) plus atoms [1 − min(α,β)]δ₀ + [max(α+β−1,0)]δ₁, r± = α+β−2αβ ± √(4αβ(1−α)(1−β)).
- Bulk: sine-kernel universality with O(n^{−1}) corrections; soft edge: Airy kernel on scale n^{−2/3}; hard edge:
  Bessel kernel on scale n^{−2}.
- **Large deviations**: for any compact K disjoint from [r,s], P(eigenvalue in K) < e^{−Cn}.

## Relationship to our project
- Exact statement behind the "Haar null model" of the decoder: with P_B Haar-random of rank d, O_T is *exactly*
  a Jacobi ensemble at finite D. The Weingarten moments and W₁ of the current draft are properties of this
  ensemble.
- The large-deviation bound justifies the claim that a Haar model has essentially no sub-edge eigenvalues
  (the audit's "finite-sample" caveat still applies at the Airy scale n^{−2/3}).
- A quantitative soft-tail test should compare the SYK sub-edge population with Tracy–Widom/Airy edge
  fluctuations of the matching Jacobi ensemble, not just with "zero".
