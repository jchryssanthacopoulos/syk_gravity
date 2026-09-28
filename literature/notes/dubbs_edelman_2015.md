# Dubbs & Edelman (2015) — Infinite random matrix theory, tridiagonal bordered Toeplitz matrices, and the moment problem

**Citation.** A. Dubbs, A. Edelman, *Infinite Random Matrix Theory, Tridiagonal Bordered Toeplitz Matrices, and
the Moment Problem*, arXiv:1502.04931 [math.PR]. File: `papers/dubbs_edelman_2015.pdf`.

## Main results
- The four "big" laws (Wigner, Marchenko–Pastur, Kesten–McKay, Wachter) all have Jacobi (three-term recurrence)
  parameters that are **Toeplitz with a length-1 boundary** (Table 2): Wachter has α₀ ≠ α₁ and β₀ ≠ β₁, with
  α_n, β_n constant for n ≥ 1.
- Tables 1–6: densities, Cauchy transforms, moments (Wachter moments via Narayana polynomials, Table 4), R/S
  transforms, free cumulants. Remarks: Narayana polynomials are the moments of MP and the free cumulants of Wachter.
- **Parametrization:** Wachter written with a, b ≥ 1 (MANOVA/Wishart ratio parameters), *different* from the
  (α,β) ∈ (0,1) projection-trace parametrization (conversion in Kunisky App. A.1).
- Proposed algorithm for the finite moment problem: Lanczos on the given moments, close the continued fraction
  with a Toeplitz tail. The authors note it "may have issues … in that atoms may emerge."

## Relationship to our project
- Cited by the current draft for its "bordered-Toeplitz" decoder density (its eqs. borderedtoeplitz–P4). The
  statement that free Wachter is length-1 bordered Toeplitz is supported here.
- **Limitation for our use:** a bordered-Toeplitz law with finite boundary length has absolutely continuous
  support exactly [α∞ − 2β∞, α∞ + 2β∞] = [λ₋, λ₊] plus at most finitely many atoms. It cannot represent a
  *continuous* sub-edge population, so it cannot be the "closed form" of the soft tail as claimed.
- The "atoms may emerge" warning applies to the draft's P₄ roots.
