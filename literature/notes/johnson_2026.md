# Johnson (2026) — Fortuitous Chaos, BPS Black Holes, and Random Matrices

**Citation.** C. V. Johnson, *Fortuitous Chaos, BPS Black Holes, and Random Matrices*, arXiv:2601.17122.
File: `papers/johnson_2026.pdf`.

## Main question
Is there a universal random-matrix model for the BPS sector (and nearby states) of extended JT supergravity,
capturing "fortuitous BPS chaos"?

## Physical system
Double-scaled Wishart-type matrix models H = M†M, M of size (N+Γ)×N, with Γ exact zero modes (BPS states).

## Key equations
- Leading density ρ₀(E) = Γδ(E) + (Γ̃/(2πħ)) √(E−E₀)/E (eq. 1): BPS delta plus a continuum starting at gap E₀,
  peaked at 2E₀, falling as 1/√E.
- Effective coupling ħ̃ = ħ√E₀/Γ̃ (eq. 2): the BPS sector is a Γ×Γ matrix model in its own right.
- Genus corrections R_g(E), W_{g,1}(z) (eqs. 4–6, 13): interpolate between Bessel (E₀⁰) and Airy (highest E₀
  power) intersection numbers; spectral curve x̄ = z² − E₀, y = z/(2(z²−E₀)) (eq. 12).
- String equation (eq. 19); universal solution u(x) = Γ̃²/x² − ħ²/(4x²) (eq. 25); exact Bessel-kernel density
  ρ(E) = (µ²/4ħ²)[J_Γ² + J_{Γ+1}² − (2Γ/ξ)J_ΓJ_{Γ+1}] (eq. 26).
- Numerical check: Marchenko–Pastur law (eq. 36) with λ± = Γ̃/2 + 1 ± √(Γ̃+1) and double-scaled edge density
  (eq. 37); Figs. 4–5.

## Results / remarks relevant here
- Eigenvalues can tunnel into the classical gap (non-perturbative), mentioned at the end.
- The non-BPS "cloud" above E₀ is interpreted as states "on the verge" of becoming BPS (chaos invasion).

## Code in this repo
`src/reproduce_johnson_fig4.py`, `src/reproduce_johnson_fig5.py` reproduce Figs. 4–5 (outputs in
`results/figures/johnson_fig*.{pdf,png}`). Conventions checked against eqs. (26), (36), (37).

## Relationship to our project
- The "energy-side" model the current draft contrasts with the decoder ("one hard edge vs two").
- The current draft's claim that the decoder's soft tail is the image of ρ_H(E) near the gap relies on this
  kind of near-BPS density. The audit and our own check indicate tail states are instead nearly *harmonic*
  (almost exactly in the old BPS space), so this link is currently unsupported.
- Johnson's model says nothing about overlaps/eigenvectors; it cannot by itself predict decoder spectra.
