# Bouchard (2024) — Les Houches lecture notes on topological recursion

**Citation.** V. Bouchard, *Les Houches Lecture Notes on Topological Recursion*, arXiv:2409.06657 [math-ph].
File: `papers/bouchard_2024.pdf`. (Read: Secs. 1 and 3.1–3.3; Airy-structure sections skimmed.)

## Content used here
- Spectral curve S = (Σ, x, ω_{0,1} = y dx, ω_{0,2} = Bergman kernel) (Def. 3.1); admissibility (Def. 3.2);
  Airy (x = z²/2, y = z) and Bessel (x = z²/2, y = 1/z) curves (Examples 3.3–3.4).
- Topological recursion (Def. 3.16, eqs. 3.23–3.24); recursion kernel K = ∫ω_{0,2}/(ω_{0,1}(z) − ω_{0,1}(σ(z)))
  (Remark 3.17, eq. 3.25).
- **For 2g−2+n > 0, ω_{g,n} has poles only at ramification points, of order at most 6g−4+2n** (after Def. 3.16).
- Connections to Hermitian matrix models via loop equations (Sec. 4) and to enumerative geometry (Secs. 3, 5).

## Relationship to our project
- Background for the current draft's Appendix A (W₁ of the Wachter curve). The pole-order check there
  (order 4 for ω_{1,1}) is consistent with the bound here.
- Normalization: TR correlators are defined relative to a chosen matrix size/ħ. The draft "fitted" an overall
  −1/a² factor. This is plausibly the ratio between expanding in 1/d² (Jacobi matrix size d) and in 1/D²
  (d = aD), up to sign conventions. It should be derived, not fitted.
- W₁ itself has been independently verified against exact Weingarten moments for k ≤ 5 (this audit) and k ≤ 6
  (earlier audit). That establishes it for the **Haar/Jacobi null model only**.
