# Chen, Colin-Ellerin, Mamroud & Papadodimas (2026) — Chaos of Berry curvature for BPS microstates

**Citation.** Y. Chen, S. Colin-Ellerin, O. Mamroud, K. Papadodimas, *Chaos of Berry curvature for BPS
microstates*, arXiv:2604.23287 (CERN-TH-2026-072). File: `papers/chen_colin-ellerin_mamroud_papadodimas_2026.pdf`.
(Reading depth: Secs. 1, 5, 6 read closely; Secs. 3, 4, 7 skimmed — 138 pp.)

## Main question
Is the non-Abelian Berry curvature of the BPS subspace, under deformations of the couplings, a random matrix
for black-hole microstates but not for horizonless/monotone states?

## Physical system
D1/D5 CFT and N=4 SYM (monotone sectors; curvature non-random, often zero); one-flavor N=2 SYK
(fortuitous; numerics up to N=18); N=2 super-JT gravity (analytic, Δ→∞).

## Key definitions / equations
- Berry curvature F_{μν,ab} = Σ_{m∉H_n} i⟨b|∂_[μ H|m⟩⟨m|∂_ν] H|a⟩/(E_n−E_m)² (eq. 1.1).
- Conjugation deformations Q → M Q M^{−1} preserve the BPS count; for commuting simple Λ_μ, Λ_ν,
  F_{μν} = i[Λ̂_μ, Λ̂_ν] with Λ̂ = P_BPS Λ P_BPS (eqs. 2.25, 5.18) — an LMRS-type observable.
- N=2 SYK: Q = Σ C_ijk ψ_iψ_jψ_k (eq. 5.2), R = N_ψ (eq. 5.4), particle–hole P (eq. 5.5).
  **Table 2** BPS degeneracies: N even: 2·3^{N/2−1} at r = N/2, 3^{N/2−1} at r = N/2 ± 1; N odd:
  3^{(N−1)/2} at r = (N±1)/2 (plus O(1) at (N±3)/2 when N ≡ 1 mod 4).
- Super-JT: n_BPS(R) = cos(πR) e^{S0} (eq. 6.4); Tr Λ̂² = C_S^{−2Δ} e^{S0} (eq. 6.8).

## Main results
1. SYK (Sec. 5.1–5.2): Berry curvature and i[Λ̂_μ,Λ̂_ν] have GUE nearest-neighbour statistics (N=18),
   GUE-like number variance (N=16), O(1)-looking Thouless time (N ≤ 16). Fig. 5a: Λ̂ for
   Λ = ψ₁ψ̄₂ + ψ₂ψ̄₁ at N=18 still shows **imprints of the UV spectrum** {0, ±2} of Λ; the Thouless time of Λ̂
   is longer than that of the commutator.
2. Super-JT (Sec. 6.1): at Δ→∞ crossings are suppressed by e^{−4Δ log(1+√2)} (p. 60, citing LMRS);
   non-crossing diagrams factorize ⇒ Catalan numbers / semicircle (eqs. 6.11–6.14); cylinder = annular
   non-crossing pairings (eqs. 6.16–6.19); connected SFF (eq. 6.23) → 2t/π. At Δ→0 crossed and uncrossed
   diagrams have equal weight ⇒ Gaussian (2k−1)!! (p. 60; footnotes 70, 72).
3. Sec. 6.2: at Δ→∞, Λ̂_μ, Λ̂_ν are asymptotically free and F_{μν} has the **free commutator of two semicircles**
   spectrum (eqs. 6.26–6.34); ramp with O(1) Thouless time. Footnote 72: for light operators F vanishes and
   there is no LMRS chaos.
4. Footnote 69: in double-scaled SYK each chord crossing is suppressed by a constant factor, a toy model for
   finite Δ.
5. Sec. 7: Chern numbers of the BPS bundle over SYK moduli space grow exponentially with N (eq. 7.25; skimmed).

## Relationship to our project
- The current draft's description of this paper (free commutator of two semicircles, chord technology) is
  accurate. Its further claim that "our role of Δ is played by D = e^{S0}" is **not** supported here:
  crossings on the disk are controlled by Δ.
- Two different "free" limits must be kept apart: (i) LMRS/Berry freeness from Δ→∞ (semicircular Ô);
  (ii) Haar freeness of the BPS subspace relative to a fixed projector (Wachter). At a→0 the Wachter law
  (recentred, rescaled by √a) also becomes a semicircle, but with width ∝ √a rather than ∝ S0^{−Δ}.
- The "UV imprint" in Fig. 5a is a close analogue of the wall enhancement (pile-up at λ=0,1) in our decoder
  spectrum, whose UV operator Π_T has spectrum {0,1}.
- Table 2 independently confirms the BPS dimensions used in our code (e.g. N′=13, P=6 → 729).
