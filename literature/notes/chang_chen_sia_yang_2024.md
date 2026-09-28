# Chang, Chen, Sia & Yang (2024) — Fortuity in SYK models

**Citation.** C.-M. Chang, Y. Chen, B. S. Sia, Z. Yang, *Fortuity in SYK models*, JHEP 08 (2025) 003,
arXiv:2412.06902 (v3). File: `papers/chang_chen_sia_yang_2024.pdf`.
(Reading depth: Secs. 1, 3.2, 3.4, 6 read; Secs. 2, 4, 5 skimmed.)

> **Citation warning.** The current draft cites arXiv:2412.06902 as "C.-M. Chang, Y.-H. Lin, *Fortuity in the
> N=2 SYK model*". That is wrong: the authors are Chang, Chen, Sia, Yang and the title is *Fortuity in SYK models*.

## Main question
How does fortuity (Chang–Lin) manifest in supersymmetric SYK, and what makes a supercharge "generic"?

## Physical system
N=2 SYK, Q = Σ C_{i1..iq} ψ_{i1}…ψ_{iq} (eq. 1.2) as wedge product with a q-form on Λ^•(C^N) (eq. 1.7);
two-flavor models (Heydeman–Turiaci–Zhao) and symmetrized variants with monotone states.

## Key results
1. **R-charge concentration**: for generic C, rank M_p(C) = D̃_p for p < (N−q)/2 (eqs. 1.8–1.9), so all BPS
   states lie in N/2 − q/2 ≤ p ≤ N/2 + q/2; all are fortuitous.
2. **Supercharge chaos conjecture** (Conjecture 1): a generic q-local supercharge has BPS states concentrated in
   one R-charge per irreducible complex and is approximated near them by the Turiaci–Witten N=2 ensemble.
3. Symmetrized two-flavor model Q = Σ C_ijk ψ_iψ_(j χ̄_k) (eq. 3.18) commutes with V = Σψ_iχ_i (eq. 3.19);
   V^k|Ω⟩ is a monotone BPS tower forming a spin-N/2 su(2) representation (eqs. 3.20–3.23).
4. LMRS test (K=2 two-flavor, N=4, N_ψ=N_χ=4; 1820 BPS states, 15 monotone): with
   O = ψ¹₁ψ̄²₁ + h.c. (eq. 3.50), Ô_m is evenly spaced/degenerate, Ô_f has GUE spacings and a long ramp
   (Figs. 6–7). Information entropy of typical fortuitous states ≈ random-state value (eq. 3.54);
   entanglement entropy larger for fortuitous than monotone states (eq. 3.56).
5. "Following N" / chaos invasion in SYK (Sec. 4); sparse couplings break concentration only when very sparse (Sec. 5.1).
6. Discussion: N^{−q} corrections could distinguish fortuitous states from Haar-random ones; "fragility of
   monotonous states" = sensitivity of monotone states to coupling perturbations.

## Relationship to our project
- Source of the one-flavor model and its BPS counting (D̃_p) used throughout our code.
- Source of the two-flavor monotone tower used in CV Sec. 6.2 and (confusingly) in the current draft's
  "chaotic vs deterministic" paragraph. The current draft mixes one-flavor and two-flavor statements.
- **Terminology clash:** "fragility" here means sensitivity of monotone states to perturbations of C; in CV it
  means *metric* fragility of uplift. Our documents should say "metric fragility" explicitly.
- The suggestion that fortuitous states differ from Haar-random states at subleading order is directly
  relevant to the non-Wachter deviations of the decoder.
