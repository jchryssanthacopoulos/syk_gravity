# Fu, Gaiotto, Maldacena & Sachdev (2016) — Supersymmetric SYK models

**Citation.** W. Fu, D. Gaiotto, J. Maldacena, S. Sachdev, *Supersymmetric Sachdev-Ye-Kitaev models*,
Phys. Rev. D 95, 026009 (2017) [Addendum: 95, 069904], arXiv:1610.08917 (v3). Journal ref from CV's
bibliography. File: `papers/fu_gaiotto_maldacena_sachdev_2016.pdf`. (Read: Secs. I, II (summary), V; N=1
numerics skimmed.)

## Main content
- N=1 model: Q = iΣC_ijkψ^iψ^jψ^k with Majorana fermions (eq. 1.1), H = Q² (eq. 1.4). At large N, unbroken SUSY
  with Δ_f = 1/6 and a bosonic partner Δ_b = 2/3. In the exact theory SUSY is broken non-perturbatively
  (E₀ ~ e^{−αN}). Generalization to q̂-body supercharges: Δ_ψ = 1/(2q̂).
- **N=2 model (Sec. V):** complex fermions, Q = iΣC_ijkψ^iψ^jψ^k, Q̄ = its conjugate (eq. 5.1), Q² = Q̄² = 0,
  H = {Q, Q̄} (eq. 5.3), ⟨C C̄⟩ = 2J/N² (eq. 5.4). Fermions have R-charge 1/q̂ and dimension Δ_ψ = 1/(2q̂)
  (chiral primary: R = 2Δ).
- **Z_q̂-twisted Witten index** W_r = Tr[(−1)^F e^{2πirQ_R}] = e^{iNπ(r/q̂ − ½)}(2 sin(πr/q̂))^N (eq. 5.5).
  Maximal at r = (q̂±1)/2, with log|W| = N log(2cos(π/(2q̂))) (eq. 5.6), which equals the large-N ground-state
  entropy. For q̂=3 the entropy density is ½ln3 per complex fermion.
- **Exact ground-state degeneracies** for q̂=3 (eq. 5.7): N even: D(N,0) = 2·3^{N/2−1}, D(N,±1/3) = 3^{N/2−1};
  N ≡ 3 mod 4: D(N,±1/6) = 3^{(N−1)/2}; N ≡ 1 mod 4: D(N,±1/6) = 3^{(N−1)/2}, D(N,±1/2) = 1 or 3. Zero otherwise.
- N=2 super-reparametrization and super-Schwarzian (Sec. V.A onward).

## Relationship to our project
- **Primary source** for the one-flavor N=2 SYK model, the index, and the BPS dimensions used throughout (CV, CCSY,
  Berry Table 2 all quote eq. 5.7). This confirms d(12→13) = 3^6 = 729 etc. It also shows the draft's attribution of
  the ½ln3 entropy density to a "Lagrange–Bürmann proof" is unnecessary: it is FGMS eq. (5.6).
- Supplies Δ_ψ = 1/(2q̂) ⇒ Δ(n_i − ⟨n_i⟩) = 1/q̂ = 1/3 for the decoder's UV operator (research Q1).
- Coupling normalization ⟨|C|²⟩ = 2J/N² differs from CV's J(q−1)!/N^{q−1} and from the unit variance used in most
  repo code. This is irrelevant for BPS subspaces (scale invariance) but matters for any energy comparison
  (e.g. E_b in the soft-tail identity).
