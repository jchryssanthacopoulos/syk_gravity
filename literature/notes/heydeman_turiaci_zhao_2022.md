# Heydeman, Turiaci & Zhao (2022) — Phases of N=2 Sachdev–Ye–Kitaev models

**Citation.** M. Heydeman, G. J. Turiaci, W. Zhao, *Phases of N = 2 Sachdev-Ye-Kitaev models*, arXiv:2206.14900
(v2, Oct 2022). File: `papers/heydeman_turiaci_zhao_2022.pdf`. (Read: Sec. 1 and model definitions; the
Schwinger–Dyson analysis in Secs. 2–3 skimmed.)

## Main results
1. **FGMS N=2 SYK at non-zero background R-charge** (Sec. 2): the IR has an emergent super-reparametrization
   symmetry broken to SU(1,1|1). The fermion scaling dimension Δ **depends on the background charge** (eq. 2.27;
   Δ = 1/(2q̂) only at zero charge), and the two-point coefficient is not fixed by the IR at nonzero charge.
2. Critical charge beyond which Δ violates unitarity: high→low entropy phase transition to a gapped free-fermion
   phase (Sec. 2.2). Interpreted as the end of black-hole stability (Schwinger pair production / BF-bound analogy).
3. Luttinger–Ward relation (eq. 2.50), grand potential (eq. 2.90), bilinear operator spectrum (Sec. 2.3).
4. **New multi-flavour model** (Sec. 3): Q = iC_ijk ψ^iψ^jχ^k (eq. 3.1), with U(1)_ψ × U(1)_χ and a supersymmetric
   flavour symmetry Q_F = Q_ψ − 2Q_χ (eq. 3.2). Integer charges and a non-vanishing index. This is the
   "two-flavor model" that CCSY and CV later symmetrize to obtain monotone states.

## Relationship to our project
- Origin of the two-flavor model used as a control in CV §6.2 and CCSY §3.
- **Caveat for research Q1:** LMRS-type predictions use the conformal dimension of the fermion bilinear. In our
  decoder the BPS sector sits at fixed fermion number near half filling (small R-charge, |Q_R| < 1/2, FGMS eq. 5.7; zero background charge), where
  Δ_ψ = 1/(2q̂) should hold, so Δ(n − ⟨n⟩) = 1/q̂ = 1/3. If one moves off the central sector, eq. (2.27)'s charge
  dependence must be included.
