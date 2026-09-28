# Penington, Shenker, Stanford & Yang (2019) — Replica wormholes and the black hole interior

**Citation.** G. Penington, S. H. Shenker, D. Stanford, Z. Yang, *Replica wormholes and the black hole interior*,
arXiv:1911.11977 (v2, Apr 2020); JHEP 03 (2022) 205 as cited in the project draft (journal ref not verified from
the PDF). File: `papers/penington_shenker_stanford_yang_2019.pdf`.
(Read: Secs. 1–2 closely and App. D; Secs. 3–7 skimmed.)

## Main question
Derive the island rule / Page curve from the gravitational replica trick, and explain how "orthogonal" bulk
states acquire small overlaps.

## "West-coast" model (Sec. 2)
- JT gravity plus an **EOW brane** with tension µ (eqs. 2.1–2.4) carrying k orthogonal internal states entangled
  with a reference R: |Ψ⟩ = k^{−1/2} Σ_i |ψ_i⟩_B|i⟩_R (eq. 2.5).
- Naive overlaps ⟨ψ_i|ψ_j⟩ ∝ δ_ij (eq. 2.10); replica wormholes give |⟨ψ_i|ψ_j⟩|² ∝ δ_ij + e^{−S0} (eqs. 2.13–2.16),
  interpreted as ⟨ψ_i|ψ_j⟩ = δ_ij + e^{−S0/2}R_ij with random R (eq. 2.17).
- Topology weighted by e^{S0χ}; Z_n ∝ e^{S0} for connected n-boundary geometries (eq. 2.14).
- **Planar resolvent** Schwinger–Dyson equation (eqs. 2.22–2.28). In JT, Z_n = e^{S0}∫ds ρ(s) y(s)^n
  (eq. 2.32), leading to λR = k + ∫ds ρ(s) w(s)R/(k − w(s)R) (eq. 2.34). In a microcanonical window this gives
  a Marchenko–Pastur entanglement spectrum and S(R) = min(log k, S_BH) (eq. 2.7).
- Footnote 15: the analysis was inspired by free-probability results.
- App. D: an explicit ensemble dual in which the brane states are random vectors with i.i.d. Gaussian coefficients
  C_{i,a} (eq. D.1), coupled to the JT matrix integral.

## Relationship to our project
- Cited by the draft as the template for a "two-brane-species" model of the decoder. PSSY has **one** species of
  EOW brane with k flavours and derives Marchenko–Pastur statistics of a Gram matrix. Nothing in PSSY supports a
  second species labelled by Fock occupation patterns (the draft's "t-branes"), or identifying e^{S0} with the
  ambient Hilbert-space dimension.
- The draft's claim "PSSY's Marchenko–Pastur law is the b→0 limit of Wachter" is at best an analogy, and needs an
  eigenvalue rescaling.
- The BPS analogue done properly is Boruch–Iliesiu–Yan (see that note).
