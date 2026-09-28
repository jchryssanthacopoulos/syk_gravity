# Chryssanthacopoulos & Vegh (2026) — Fortuity and fragility in supersymmetric SYK

**Citation.** J. Chryssanthacopoulos, D. Vegh, *Fortuity and fragility in supersymmetric SYK*,
arXiv:2608.12160 [hep-th] (QMUL-PH-26-29), 12 Aug 2026. File: `papers/chryssanthacopoulos_vegh_2026.pdf`.

**Role in this project.** This is the project's own published paper and the *direct parent* of the
current draft (`research/pdfs/gravitational_dual_fragility.pdf`). It defines the decoder, the
fidelity/dressing-cost diagnostics, and the operator O_T whose spectrum the current draft studies.
Its Discussion ends by leaving "a gravitational construction of the decoder" to future work — that
is the stated goal of the current draft.

## Main question
Can BPS states of N=2 SYK be continued ("uplifted") from N to N+1 fermions, and what does the
Hilbert-space *metric* (not just cohomology) say about how well they survive? Does this distinguish
fortuitous (chaotic) from monotone (protected) BPS sectors?

## Physical system
- One-flavor N=2 SYK: N complex fermions, H_N = Λ^•(C^N), supercharge Q_N = Σ C_{i1..iq} ψ_{i1}…ψ_{iq}
  (eq. 4), generic complex Gaussian couplings, q odd (mostly q=3). BPS = ker Q ∩ ker Q† ≅ H^p(Q_N) (eq. 6),
  Hodge decomposition eq. (7).
- Solvable controls: Chen's single-matrix model Q = Tr Ψ³ (Sec. 6.1) and the symmetrized two-flavor
  model Q = Σ C_ijk ψ_i ψ_j χ†_k with C symmetric in (j,k) (eq. 82), which has a monotone tower V^n|Ω⟩,
  V = Σ ψ_i χ_i, [Q,V]=0 (eqs. 83–84).

## Key definitions
- Z_q-twisted index χ(r) (eqs. 9–11): for q=3, χ(r) = (2/3)·3^{N/2}(−1)^r cos(π(N−2r)/6).
- Fortuity window (N−q)/2 ≤ p ≤ (N+q)/2 (eq. 12); level ℓ = 2p − N (eq. 29); exceptional boundary
  states at ℓ = ±q (eqs. 14–16).
- Adding one fermion: Ψ = α + ψ_{N+1} β (eq. 17); Q_{N+1} = Q_N + C′ψ_{N+1}, C′ = Σ C_{ij,N+1}ψ_iψ_j
  (eq. 19); block form (eq. 21); closedness ⇔ Q_N α = 0, C′α = Q_N β (eq. 22).
- Long exact sequence (eq. 24), connecting map δ[α] = [C′α] (eq. 25); dimension formula (eq. 26) and,
  when connecting maps vanish, dim H^p(Q_{N+1}) = dim H^p(Q_N) + dim H^{p−1}(Q_N) (eq. 27).
- Cohomological walks / binary strings (eq. 28); two uplift channels (eqs. 30–31).
- **Decoder** D_T := Π_T P_{B^P_{N′}} : B^P_{N′} → H^p_N (eq. 42): read off the T-occupation component of an
  enlarged BPS state. Exact-lift locus A^p_{N;T} = B^p_N ∩ im D_T (eq. 44).
- K_T := D_T D_T† on H^p_N (eq. 45); SVD (eqs. 46–47); eigenvalues λ_j = σ_j² ∈ (0,1].
- **State-dependent** decoder spectral measure dμ_{α,T} for an old BPS state α (eq. 49).
  Uplift fidelity F_T(α) (eq. 50), bare fraction f_{0,T}(α) (eqs. 54–56), inverse moment m_{−1,T}(α)
  (eq. 58–60), dressing cost κ_T(α) = m_{−1} − 1 (eq. 61), projected dressing cost (eq. 62).
- O_T := D_T† D_T = P_B Π_T† Π_T P_B on the enlarged BPS space (eq. 101); for one added mode
  O↓ = P_B(1−n_{N′})P_B, O↑ = P_B n_{N′} P_B, O↓ + O↑ = 1_B (eq. 102).

## Main results
1. Inside the fortuity window at least one of the two connecting maps vanishes (eqs. 32–35), so every
   nonzero BPS class has a one-step uplift; iterating, every class can be uplifted to arbitrarily large N
   along a walk (Sec. 4.2). Uplift is non-canonical (coset ambiguity eq. 40; harmonicity not preserved,
   eqs. 36–39).
2. Single-matrix model: bare staircase embedding E_N(α) = α ⊗ T_row is isometric and BPS-preserving
   (eqs. 73–78), so F=1, f_0=1, κ=0 (eq. 80).
3. Two-flavor tower: V^n_{N+1}|Ω⟩ = |V^n_N⟩⊗|00⟩ + n|V^{n−1}_N⟩⊗|11⟩ (eq. 85); F_00 = F_11 = 1;
   κ_00 ≤ n/(N+1−n) (eq. 87), κ_11 ≤ (N−n)/(n+1) (eq. 88). Fortuitous states: minimum fidelity U_T (eq. 89,
   Fig. 6), dressing costs grow with N (Figs. 8–11).
4. One-flavor model: generic-position counting inside the *constrained rooms* R_↓ = ker Q_N, R_↑ = ker Q_N†
   (eqs. 93–95) reproduces the exact-lift dimensions (Table 1) except at (9,3). Rank bound (eq. 98) and
   Conjecture A.1 give codim > 0 for odd N ≥ 15 and A = {0} for odd N ≥ 19 (eq. 100) — *conditional*
   (Conj. A.1, vanishing connecting maps, generic position). A direction is already lost at N = 11.
5. **BPS chaos of the decoder**: for N′=12, P=6, 500 realizations, the interior spectrum of O↓ has
   ⟨r⟩ = 0.599 ± 0.001 (eq. 105), consistent with GUE (benchmarks eq. 104; GUE surmise 0.60266).
6. Appendix A: ordinary BPS multiplicities = strip-confined ±1 walks, h^p_o = D̃_p − D̃_{N−p−q}
   (eqs. 125, 130), with D̃_p = Σ_n (−1)^n C(N, p−nq) (eq. 114) and generating function (1+x)^N/(1+x^q)
   (eq. 115). Conditional on Conjecture A.1 (rank formula, eq. 120), verified numerically q = 3…9, N ≤ 14.

## Numerical methods
Exact diagonalization / nullspace computations of the wedge operator on Λ^•(C^N); up to 500 disorder
realizations for ⟨r⟩; five realizations for U_T; single realizations for κ plots. Code in this repo:
`src/reproduce_table1.py` (Table 1), `src/decoder_chaos.py` (Fig. 13 / eq. 105).

## Limitations (as stated or evident)
- Asymptotic obstruction (eq. 100) rests on Conjecture A.1 and generic position.
- κ-data are single-realization and basis-dependent (state-resolved values depend on the chosen harmonic basis).
- ⟨r⟩ only at one (N′,P); Thouless time and N-scaling not studied.
- The two-flavor fortuitous analysis reaches only N = 8.

## Relationship to our project
- Defines everything the current draft uses. **Caution:** the current draft redefines
  "uplift fidelity" as F_T = (1/d) Tr O_T (uniform average over the *enlarged* BPS space) whereas here
  F_T(α) is a *state-dependent* quantity for an *old* BPS state α (eq. 50). These are different objects.
- Eq. (27) implies d = dim B^p_N + dim B^{p−1}_N for the enlarged sector, i.e. a canonical *reference
  subspace* B⁰ = B^p_N ⊕ (B^{p−1}_N ∧ e_{N′}) of the same dimension d on which O↓ has spectrum {0,1}.
  This is the natural starting point for understanding the non-Wachter structure of the decoder
  (see `docs/research_questions.md`, Q2).
- Table 1's "generic position inside constrained rooms" logic is the right tool for the exact λ=0,1 atoms
  (the current draft instead calls them a "finite-N resonance").
- The paper explicitly does **not** claim a gravity dual of the decoder.
