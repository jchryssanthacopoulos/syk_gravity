# Boruch, Iliesiu & Yan (2023) — Constructing all BPS black hole microstates from the gravitational path integral

**Citation.** J. Boruch, L. V. Iliesiu, C. Yan, *Constructing all BPS black hole microstates from the gravitational
path integral*, arXiv:2307.13051 (v2, Jul 2024). File: `papers/boruch_iliesiu_yan_2023.pdf`.
(Read: Secs. 1–3.5 closely; 3.6–4 and appendices skimmed.)

## Main question
Can the gravitational path integral prepare an explicit basis of BPS black-hole microstates, and does the state
count reproduce the Gibbons–Hawking degeneracy, including non-perturbative corrections?

## Physical system
N=2 JT supergravity (N=2 super-Schwarzian) with probe matter, in the zero-temperature (BPS) limit. Motivated
by 1/16-BPS AdS₅ black holes and 4d N=8 BPS black holes.

## Key definitions / equations
- BPS Hartle–Hawking wavefunction Ψ_j(ℓ, a) in super-Liouville variables (eqs. 3.7–3.9); Z_j = e^{S0}cos(πj)
  (eq. 3.15).
- Prepared states |q_i⟩ = O_i|HH⟩ with infinite Euclidean evolution on both sides of a (not necessarily BPS)
  matter insertion (eqs. 2.4–2.6). Equivalently K_i = lim e^{−βH}O_i e^{−βH}, an LMRS operator (eq. 2.5).
- Disk two-point function ⟨2pt⟩_{j,disk} (eq. 3.19); cylinder with a matter propagator (eqs. 3.25–3.29);
  n-boundary "pinwheel" wormholes Z_{j,O_n} ∝ e^{S0(2−n)}cos²(πj)(…)^n (eqs. 3.30–3.32, 3.43).
- Gram matrix M_ij = ⟨q_i|q_j⟩ (eq. 3.37); resolvent Schwinger–Dyson equation over planar pinwheels
  (eqs. 3.42–3.46).

## Main results
1. The resolvent gives a **Marchenko–Pastur** eigenvalue density for M with edges λ± = y(√K ± e^{S_j})²
   (eq. 3.47), e^{S_j} ≡ e^{S0}cos(πj). The rank of ρ = Σ|q_i⟩⟨q_i| saturates at e^{2S_j} = (dim H_BPS)²
   (eqs. 2.10, 3.50–3.51). Null states appear for K > e^{2S_j}.
2. **Only genus-zero wormholes contribute** (Sec. 3.4). Geometries containing closed geodesics that avoid all
   matter propagators (trumpets, higher genus) vanish at β→∞. The only extra corrections come from single
   supersymmetric defects, which shift the dimension to (Z_BH^{β→∞})² (orbifold contributions). The sum over
   geometries is convergent.
3. **Haar-random-state dual** (Sec. 3.5, eqs. 3.68–3.72): |q_i⟩ ∝ Σ_a C_ia |BPS_a⊗²⟩ with i.i.d. complex Gaussian
   C_ia reproduces all nonzero gravitational contributions, including non-planar ones. The vanishing of trumpets
   is what makes this simple Gaussian model exact (fn. 15; contrast PSSY, where the random states must be
   coupled to the JT matrix integral).
4. The rank has vanishing standard deviation (Sec. 3.6): the BPS count is exact, not ensemble-averaged.
5. Any state in H_BPS ⊗ H_BPS can be reconstructed from the |q_i⟩ (Sec. 4).

## Limitations (stated)
Non-planar geometries and full matter loop effects neglected in the resolvent (Sec. 3.7). Works in N=2 JT as a
proxy for N=4 JT.

## Relationship to our project
- This is the correct, worked-out version of what the current draft gestured at ("PSSY sum over EOW branes in the
  E=0 sector of N=2 super-JT"). It **contradicts** three of the draft's ingredients:
  1. the dimension parameter is **e^{S_j} = e^{S0}cos(πj) = BPS degeneracy** (our d), not the Fock-sector D;
  2. in the BPS sector **higher-genus / trumpet contributions vanish**, so there is no gravitational "handle
     tower ∝ D^{−2g}" of the type the draft attributes to gravity;
  3. states are prepared by **matter-operator insertions** with Δ-dependent weights (eq. 3.43), not by a second
     brane species labelled by slot occupation patterns.
- It supports a *different* statement that is useful: gravity predicts that BPS states probed by simple operators
  look like **Haar-random vectors inside the BPS space** (eq. 3.68). For the decoder, this suggests modelling the
  matrix elements ⟨ψ_i|Π_T|ψ_j⟩ (i,j ∈ BPS) as ETH-like random matrices with variance set by an LMRS two-point
  function. The draft's model is instead Haar randomness of the BPS subspace inside the ambient Fock space.
- Relevant to research Q1 and Q4 (`docs/research_questions.md`).
- Related follow-up to acquire: Boruch, Lin, Yan, *Exploring supersymmetric wormholes in N=2 SYK with chords*,
  JHEP 12 (2023) 151 (cited in Berkooz–Mamroud).
