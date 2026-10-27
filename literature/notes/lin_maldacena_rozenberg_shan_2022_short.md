# Lin, Maldacena, Rozenberg & Shan (2022) — Holography for people with no time

**Citation.** H. W. Lin, J. Maldacena, L. Rozenberg, J. Shan, *Holography for people with no time*,
SciPost Phys. 14, 150 (2023), arXiv:2207.00407 (v3, 12 Mar 2023; the journal reference is as given in CV 2026 ref. [26]
and G. Lin 2026 ref. [39], not on the PDF). File: `papers/lin_maldacena_rozenberg_shan_2022_2.pdf`. This is the short companion of
`lin_maldacena_rozenberg_shan_2022.md` (arXiv:2207.00408, the "long LMRS"). Read in full, 2026-10-27.

**Role in this project.** It gives the physics behind the LMRS operators Ô = P O P in compact form, plus the only N=2 SYK numerics
of LMRS that are quantitative and cheap to cite:
- the gap formula (33), with Fig. 9;
- the projected two-point functions in Table 1.

Of these, eq. (33) matters most for us. It is a citable statement that the continuum gap in an R-charge sector closes at
|j| = ½, which is our window edge and the source of the D10 resonant regime.

## Main question
In the zero-temperature limit of an extremal N=2 black hole, H = 0 on the ground states. What can the AdS₂ throat
still compute (correlators of projected operators, wormholes between two extremal black holes), and what does that say
about the microstates?

## Physical system
- N=2 JT gravity, whose boundary mode is the N=2 super-Schwarzian, with probe bulk matter.
- Compared with N=2 SYK (q̂ = 3) at N = 16 complex fermions.

## Important definitions
- Ô = P O P, with P = lim_{β→∞} e^{−βH} the projector on the zero-energy (BPS) states. Hats always mean "evolved for a
  very long Euclidean time" (text below eq. 13).
- The two-sided state is a superposition of wormhole lengths ℓ. The two-point function is a matrix element of e^{−Δℓ}
  (eq. 4).
- Units C = ½ in the action (5), with t = 2Cu (6). C is the Schwarzian coupling.
- The zero-energy partition function is Z(β = ∞) = e^{S₀} (9).

## Main assumptions
- Semiclassical throat with a probe-matter approximation.
- Topology changes are mentioned only as e^{−S₀} corrections (footnote 2).
- N even in the SYK comparison (footnote 9). Odd N needs a modification that the paper does not discuss.

## Main analytical results and equations
1. **Constant long-time correlators** (§2.2): ⟨O(u)O(0)⟩_β → const for u, β, β−u ≫ 1 (3).
   - Eq. (10): lim ⟨TFD(u′)|e^{−Δℓ}|TFD(u)⟩ = e^{S₀}⟨0|e^{−Δℓ}|0⟩ = e^{S₀} Γ(1+2Δ)/4^{2Δ}.
   - Eq. (11): the zero-energy wormhole has finite length, ⟨ℓ⟩ = 2(γ_E + log 4).
2. **Physical normalization** (13), (24): c_Δ² = ε^{2Δ}Γ(1+2Δ)/4^{2Δ} ∝ S_e^{−2Δ}.
   - The typical eigenvalue of Ŵ is of order c_Δ (25).
   - The factor S_e^{−Δ} comes from propagating a distance d̃ = log S_e into the interior (26). LMRS read this as: "all
     microstates must agree with AdS₂ up to that distance" (fuzzball constraint, Fig. 8).
3. **Disk/cylinder relation** (§2.3, eqs. 27–29). The cylinder two-point function in the probe approximation equals the
   disk result. So the average over couplings agrees with the average over states, to leading order.
4. **n-point functions** (§2.4). The time-independent correlators of Ô are built from the e^{−Δℓ} matrix elements.
   - TOC = OTOC for BPS (charged) fermions.
   - For heavy operators, crossed contractions are suppressed. P O P is then effectively a Gaussian random matrix on the
     BPS space, with a semicircle spectrum.
   - The detailed formulas are in the long paper (§5.3 there; see `lin_maldacena_rozenberg_shan_2022.md`).
5. **Wormholes between extremal black holes** (§2.6). Gaussian random BPS states get an end-of-the-world-brane dual
   (§2.7 interpretation).
6. **Gap per R-charge** (§3, eq. 33): E_gap(j) = (|j| − ½)²/(8C), with C = α_S J N and α_S = 0.00842…
   - α_S comes from the large-N Schwinger–Dyson equations (their ref. [34]).
   - The formula is attributed to [10] (Heydeman–Turiaci–Zhao in our collection).
   - It holds for R-charges j outside the BPS window |j| < ½.

## Numerical methods
- Exact diagonalization of N=2 SYK (q̂ = 3, N = 16), with a single realization for Fig. 9.
- Table 1 compares projected two-point functions with the Schwarzian prediction, using the large-N C = α_S N = 0.1347.
  Examples:
  - ψ̄_iψ_j at R = 0: predicted 0.0874, measured 0.079.
  - The neutral ψψ̄ OTOC/TOC ratio is about 76 %.

## Relevant figures
- Fig. 2: the wormhole-length picture.
- Fig. 8: the fuzzball depth d̃.
- Fig. 9: the gap against j at N = 16. The solid curve is (33). The data follow it outside the window, and the gap is
  zero inside.

## Limitations
- Short version: most derivations are deferred to arXiv:2207.00408.
- The SYK checks use one size (N = 16) and, for Fig. 9, one realization.
- The large-N C is used at N = 16 without a finite-N correction.

## Relationship to our project
- **Window edge, gravitational side.** In our convention the chord charge s = (degree − N/2)/p is the same as the FGMS/LMRS
  R-charge j (fermions carry 1/q̂; D11 setup). Eq. (33) then says the continuum threshold of the target sector of a
  dressing scales like (s_in − ½)²/(8C).
  - It is finite in our gapped regime (s_in > ½).
  - It vanishes in the resonant regime (s_in = ½).
  - [interpretation] This is the gravitational reason for the D10 three-regime classification and for the heavy
    K^{−1/2} tail at s_in = ½ (D11 layer 1). The full statement and the proposed test are in `literature/synthesis.md`
    §1, 2026-10-27.
- **Free-probability side.** At large Δ, P O P has a semicircle spectrum (item 4), i.e. it behaves like a free/Gaussian
  random matrix. Our decoder is a *projector* compressed between BPS projectors, and its Jacobi/MANOVA law is the
  analogue of this statement for a non-primary composite operator (see the long-LMRS note).
- **Disk/cylinder (item 3)** is the gravitational version of our use of disorder averages as state averages.
- **Normalization (13)/(24)** fixes how projected correlators scale with S_e. This is relevant if we ever compare
  absolute moments rather than shapes.
