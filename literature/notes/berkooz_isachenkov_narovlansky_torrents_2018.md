# Berkooz, Isachenkov, Narovlansky & Torrents (2018) — Towards a full solution of the large N double-scaled SYK model

**Citation.** M. Berkooz, M. Isachenkov, V. Narovlansky, G. Torrents, *Towards a full solution of the large N
double-scaled SYK model*, JHEP 03 (2019) 079, arXiv:1811.02584 (v2, 5 Feb 2019). File:
`papers/berkooz_isachenkov_narovlansky_torrents_2018.pdf` (80 pp). ("BINT" in project documents.) The journal
reference was checked against the bibliography of Boruch–Lin–Yan (their ref. [9]).

**Note written 2026-10-05.**
- Read closely: §§1–2, 3.3, 4, 5, 6.1, 8, 9 and App. A.1.
- Skimmed: §3.1–3.2 (combinatorial construction of the bi-local operator) and §6.2 (Lyapunov exponent).
- Not read: §7 (Hamiltonian reduction on U_q(su(1,1))) and Apps. B–E.

## Main question
Can all correlation functions of double-scaled (Majorana) SYK be computed exactly, at all energies, from chord
combinatorics? And which algebraic structure organizes them?

## Physical system
Majorana SYK, H = i^{p/2} Σ J_I ψ_I (1.1), in the double-scaled limit N, p → ∞ with λ = 2p²/N fixed and q = e^{−λ}.
Probes are random operators M_A of length p_A with independent Gaussian couplings (2.6–2.7).

## Important definitions
- **Chord partition function** m_k = Σ_{diagrams} q^{#intersections} (2.5). It follows from Poisson-distributed
  overlaps of index sets with mean p²/N (2.4, App. A.1).
- **Transfer matrix** T on the chord-number space |l⟩ (2.13). Then m_k = ⟨0|T^k|0⟩ (2.14). The energies are
  E(θ) = 2cosθ/√(1−q) (2.15), with q-Hermite eigenvectors (2.16) and measure dμ(θ) = (dθ/2π)(q, e^{±2iθ}; q)_∞
  (5.2). T is a position operator A + A⁺ of the Arik–Coon q-oscillator (2.20–2.21).
- **Matter rule (2.9–2.11).**
  - An H-chord crossing an M_A-chord gives q̃_A = e^{−2pp_A/N}.
  - An M_A-chord crossing an M_B-chord gives q̃_AB = e^{−2p_Ap_B/N}.
  - Only same-flavour M's pair (independent couplings).
- **Bi-local operator.**
  - A contracted pair of M's acts on the chord space. It commutes with T (§3.3, eq. 3.20), i.e. it conserves
    energy.

## Main assumptions
- Strict double scaling: overlaps are Poisson and there are no triple intersections (App. A.1, valid for overlap
  k ≪ p, p′ ≪ N).
- Random (Gaussian-coupled) matter. Footnote 4 says the method "can be adapted" for non-random operators, but
  this is not done.

## Main analytical results
1. **Two-point vertex (5.4):** √[(q̃²;q)_∞ / (q̃e^{i(±θ₁±θ₂)};q)_∞].
   - The square of this vertex, between two energy eigenstates, is exactly the Rogers–(q-)Mehler kernel
     p_r^{(q)}(x,y) with r = q̃ (Bożejko–Kümmerer–Speicher 1997, Thm. 1.10(2)).
2. **Uncrossed 4-pt function** ⟨M₁M₁M₂M₂⟩, eqs. (4.8)–(4.9).
   - It is a product of two 2-pt vertices with a shared energy.
   - It depends only on β₂ + β₄ (energy conservation of the bi-local).
3. **Crossed 4-pt function** ⟨M₁M₂M₁M₂⟩, eqs. (4.12)–(4.13).
   - It is the product of four vertices times an **R-matrix** (5.6), given by an ₈W₇ very-well-poised series.
   - Integral (Mellin–Barnes) form: (6.2)–(6.3).
   - As q → 1⁻ at low energy it reduces to the Schwarzian 6j-symbol of SU(1,1) (Mertens et al.).
   - It is related to a 6j-symbol of U_{q^{1/2}}(su(1,1)) (§7).
4. **Feynman rules (§5).**
   - Propagator e^{−Δτ E(θ)} (5.1), measure (5.2), vertex (5.4).
   - Every crossing of two matter lines carries R.
   - q̃_A = q^{l_A} defines the "dimension" l_A (5.3).
5. **Chaos (§6.2):** maximal Lyapunov exponent at low temperature, with corrections in temperature and λ.
6. **Finite-p structure (§8).**
   - The two-stage limit (N → ∞ at fixed p, then large p) and the double-scaled limit capture complementary
     terms. Moments are double series in λ and 1/p (8.1).
   - Chord diagrams with q-weights miss O(1/p) terms coming from non-Poisson overlap statistics (footnote 18).

## Important equations for us
(2.4), (2.11), (2.13)–(2.19), (3.20), (4.9), (4.13), (5.4)–(5.6), (8.1), App. (A.1).

## Numerical methods
None beyond plots of the analytic formulas (Lyapunov corrections in §6.2). Table 2 compares the chord and
Schwinger–Dyson expansions of the fermion 2-pt function at high temperature.

## Limitations
- Majorana, non-supersymmetric. No ground-state degeneracy and no BPS projector.
- Random matter only. A fixed operator (e.g. our single-site parity U) obeys different crossing rules (see below).
- Strict double scaling. At p = 3 the triple-intersection and hypergeometric-overlap corrections of §8 / App. A are
  not small.

## Relationship to our project
1. **The 2-pt vertex is the q-Mehler kernel. So "x per chord" ≡ "x^{N̂} on chord space" ≡ OU-correlated couplings.**
   - With Bożejko–Kümmerer–Speicher (Thm. 2.11) and Pluma–Speicher (Cor. 3.2), the double-scaled pair
     (H(ω), H(ω_x)) with OU-correlated couplings is a q-Gaussian pair ⟨h, h_x⟩ = x.
   - Hence E tr[f(H(ω)) g(H(ω_x))] = ⟨f(T)Ω| x^{N̂} |g(T)Ω⟩, which is BINT's 2-pt structure with q̃ → x.
   - This is the double-scaled form of our finite-N statement "chord number = Wiener chaos order"
     (`research/notes/chord_size_dictionary_2026-10-02.md`).
2. **BINT's crossed 4-pt function is *not* our four-arc OU function. The chord weights differ** [derived in
   project, `docs/derivations.md` D5].
   - On a circle with four arcs separated by M₁, M₂, M₁, M₂, an H-chord between *adjacent* arcs crosses one M-chord
     (weight q̃). A chord between *opposite* arcs crosses both (weight q̃₁q̃₂).
   - In T₄^OU(x) = E τ(P(ω)P(ω_x)P(ω)P(ω_x)), adjacent-arc chords weigh x but opposite-arc chords weigh 1. This is a
     ℤ₂ rule (write H = cA + sB, H_x = cA − sB with c² − s² = x), not a product rule.
   - So the R-matrix (5.6) cannot be used for T₄ by substituting q̃ = x.
   - The same ℤ₂ structure holds exactly for the decoder, where U = (−1)^{n_i} flips the sign of the B-terms (those
     containing site i).
3. **The algebraic target for R.** The decoder remainder R is a crossed-channel quantity (see D5). If it admits a
   closed form in double scaling, it should be a ℤ₂-twisted analogue of the R-matrix: the 6j-type crossing kernel
   of the ℤ₂ (fixed-operator) chord rule. Neither BINT nor BBNR computes such a kernel. This is a statement about
   what is missing, not a result.
4. **Finite-p corrections (§8).** These explain why double-scaled formulas are only good to a few-to-10 % at
   p = 3. They are consistent with our measured chord-weight deviations.
