# Boruch, Lin & Yan (2023) — Exploring supersymmetric wormholes in N = 2 SYK with chords

**Citation.** J. Boruch, H. W. Lin, C. Yan, *Exploring supersymmetric wormholes in N = 2 SYK with chords*,
arXiv:2308.16283 [hep-th] (v3, 19 Dec 2023); JHEP version with a Mathematica supplement.
File: `papers/boruch_lin_yan_2023.pdf`. ("BLY" in project documents. Not to be confused with Boruch–Iliesiu–Yan
2023, `boruch_iliesiu_yan_2023.md`.)

**Note written 2026-09-28.** Read in full: §§1–4, §5.1 in detail; §5.2–5.3, §6, App. A, C.2 skimmed;
App. B, C.1, C.3, D, E, F not read line by line.

## Main question

What does the zero-temperature two-sided state (the supersymmetric wormhole, i.e. the β = ∞ TFD) of N = 2 SYK look
like in the exactly solvable double-scaled limit? What does it compute: ground-state counts, zero-temperature
correlators, OTOCs?

## Physical system

- **Model:** one-flavor N = 2 SYK, the same as ours. N complex fermions, supercharge
  Q = Σ_I C_I Ψ_I with Ψ_I a product of p fermions, H = {Q, Q̄} (eq. 1.3).
  - Variance ⟨C_I C_I′*⟩ = 2 binom(N,p)^{−1} J² δ (eq. 1.4). This is a J² convention, different from FGMS/LMRS.
- **U(1)_R:** J = (1/2p) Σ[ψ_i, ψ̄_i] (eq. 1.5), normalized so that Q has charge 1. With n_i = ψ̄ψ,
  J = (N/2 − N_f)/p, the opposite sign to our j = (P − N′/2)/p; all results are even in j.
- **Double-scaling limit:** N, p → ∞ with λ = 2p²/N fixed, q = e^{−λ} (eq. 1.6).
  - The triple-scaling limit (β → ∞ and λ → 0) gives the N = 2 super-Schwarzian of LMRS.

## Important definitions

- **Chord rules** (eqs. 1.13–1.15): Q-chords carry arrows. A pair of chords is a "friend" (factor q^{−1/2}) or an
  "enemy" (±q^{+1/2}). Each diagram gets (−1)^{#intersections} q^{(#enemies − #friends)/2}.
  - A matter chord with p′ fermions (Δ = p′/p) crossing a Q-chord: q → q^Δ.
  - Matter–matter crossing: q → q^{Δ²}.
- **Conformal dimension** of W ∼ (p_O ψ's)(p_X ψ̄'s): Δ = (Δ_O + Δ_X)/2, with R-charge Δ_O − Δ_X (eq. 1.8).
  For our ψ̄_iψ_i: Δ = 1/p.
- **Two-sided chord Hilbert space** with states |n, ab, j⟩, a, b ∈ {O, X} (eqs. 2.4–2.5, 2.34).
  - n is the chord number, the discrete wormhole length, ℓ = 2λn.
  - The four "ab" labels (two qubits) are the fermionic superpartners of the Schwarzian length mode.
- **Bookkeeping label m** for "forgotten" friends and enemies (closed chords, §2.2). Its conjugate gives the U(1)_R
  charges, eqs. (2.26)–(2.27), (2.33). Relation to other papers: j = −2 j_R; LMRS uses j_LMRS = j/2.
- **Hartle–Hawking (supersymmetric-wormhole) state |Ψ, j⟩:** the state annihilated by all four bulk
  supercharges (eqs. 3.1–3.7). Its amplitudes α_n, β_n are q-Hermite polynomials.

## Main assumptions

- Double-scaling limit (p → ∞). Nothing is claimed for finite p. We test finite p = 3 ourselves below.
- The chord Hilbert space is built by a GNS-type construction. The bulk supercharges (eqs. 2.15, 2.32) satisfy
  an N = 4 algebra, checked by computer (Mathematica).
- Index saturation for the microscopic ground-state count (App. A).
- Matter operators are random p′-body operators (standard chord rules). Their normalization is fixed so that the
  n = 0 (infinite-temperature) two-point function is 1 (text after eq. 4.4).

## Main analytical results

1. **Transfer matrix with enhanced SUSY.** The empty-wormhole transfer matrix has **enhanced N = 4 SUSY**, with
   explicit bulk supercharges (2.15)/(2.32) and algebra (2.16)–(2.19).
2. **The HH ground state** in closed form (3.6)–(3.7), with norm (3.11):
   ⟨Ψ, j|Ψ, j⟩ = 1/[(q^{1±2j_R}; q²)_∞ (q²; q²)_∞].
3. **"Probability that the SUSY wormhole has zero length = BPS fraction"** (eq. 1.1, derived in 3.13–3.16).
   Double-scaled result (3.17)–(3.19):
   **D(j_R) = (q^{1+2j_R}; q²)(q^{1−2j_R}; q²)(q²; q²) = ϑ₄(iλ j_R, q)**.
   - It agrees with the double-scaled refined index (3.21)–(3.23), including non-perturbative terms e^{−(2k+1)²π²/4λ}
     that the super-Schwarzian misses.
   - Exact finite-N index: **(A.4)**.
   - Lesson (§6.3): "simple gravity" does not reproduce the exact count; its UV completion (the chords) is needed.
4. **Zero-temperature two-point function of a neutral operator** in the fixed-charge ground-state sector, at
   **finite λ** (4.5)–(4.8), derived in App. C.2:
   ⟨q^{2Δn}⟩_j = (1−q²)^{2Δ} Γ_{q²}(Δ+1)² Γ_{q²}(½ ± j_R + Δ) / [Γ_{q²}(½ ± j_R) Γ_{q²}(2Δ+1)].
   - Triple-scaling limit (4.9): (2λ)^{2Δ} cos(πj_R)/(2π) · ΔΓ(Δ)²Γ(Δ+½±j_R)/Γ(2Δ), i.e. LMRS eq. (85)/(66) with
     the fixed-sector normalization (footnote 18).
   - This identifies the Schwarzian coupling as **C = α_S N = (4λ)^{−1}**.
5. **Wormhole length.** ⟨ℓ⟩ = −2 log(1−q²) + ψ_{q²}(½ ± j_R) (4.11) is finite at finite λ. The renormalized length
   gives −ψ(½ ± j_R) in the Schwarzian limit (4.12). The small-λ expansion is asymptotic (4.13).
6. **Charged BPS two-point function** (4.14)–(4.16).
7. **Wormholes with one matter particle** (§5.1): supercharges (5.9)–(5.12), separate commuting H_L, H_R, and the
   neutral horizontal chord W_L W_R ∼ q^{Δ(n_O,tot + n_X,tot)} (5.18).
   - **The zero-temperature OTOC is set up but not evaluated** ("would give the zero temperature OTOC").
   - Generalization (5.19): P(zero length, W) = Tr(Π₀ W Π₀ W)/dim H_j.
8. **Superchord algebra** (§5.2) and coproduct (§5.3), built on n_O, n_X, which relate to operator size
   (eqs. 5.20–5.23).
9. **Discussion.**
   - Krylov/Gram–Schmidt bulk-to-boundary map; the supersymmetric wormhole has bounded Krylov complexity at
     finite λ (§6.1).
   - A sign problem prevents a "typical boundary length" picture (§6.2).
   - Remarks on higher-dimensional extremal black holes (§6.3).

## Important equations for us

- (1.6) λ = 2p²/N.
- (1.8) Δ = (Δ_O + Δ_X)/2.
- (3.11), (3.17), (A.4): BPS fraction.
- (4.7)/(4.8): finite-λ zero-temperature two-point function.
- (4.9), and C = α_S N = 1/(4λ).
- (5.9)–(5.12), (5.18): one-particle wormholes and OTOCs.

## Numerical methods

- No SYK numerics.
- Mathematica checks of the SUSY algebra and supplementary notebooks.
- Figures 1–4 are plots of the analytic formulas.

## Relevant figures

- **Fig. 1:** HH length distribution P_n.
- **Fig. 2:** D(j_R) versus j_R for several λ.
- **Fig. 3:** zero-temperature two-point function versus Δ, approaching the Schwarzian as λ → 0.
- **Fig. 4:** ⟨ℓ⟩ versus q.

## Limitations

- **Double-scaled throughout (p → ∞).** For fixed p, the λ → 0 limit of these formulas is not the true large-N
  theory. Its Schwarzian coupling is α_S = 1/(8p²) (= 0.0139 at p = 3), whereas LMRS quote 0.00842 for the
  q̂ = 3 model.
- **The OTOC and higher-point zero-temperature correlators are not computed.** Only the Hilbert-space machinery
  is given.
- **The chord sums are not sign-definite** (§6.2).
- **Matter is modeled as random p′-body operators.** Fixed, site-local operators such as n_i need their own
  crossing weights. In particular, two different single-site bilinears commute exactly, so their
  matter–matter crossing weight is 1, not q^{Δ²}.

## Relationship to our project (tested numerically 2026-09-28; see `docs/derivations.md` D2)

- **Same model, same observables.** BLY's Tr(Π_{0,j} W Π_{0,j} W) in a fixed-charge sector, normalized to the
  infinite-temperature value, is **exactly our r = Var(λ)/[b(1−b)]** for W = n_i − μ, with Δ = 1/3 (p = 3).
  - Their D(j) is our a = d/D.
  - Their (5.18) OTOC is our crossing weight X.
- **BPS fraction.** The finite-N index (A.4) reproduces all our exact BPS dimensions (every sector
  |j| < ½, N′ ≤ 16; test `tests/test_dssyk_n2.py`). The double-scaled D(j) (3.17), evaluated at λ = 18/N′,
  overshoots the exact a by 0–10 % (growing with N′). That is the size of the finite-p error:
  cos^N(π/2p) versus e^{−Nπ²/8p²}.
- **Decoder variance: a parameter-free prediction that works.**
  - Eq. (4.8) at λ = 18/N′, Δ = 1/3 matches the measured r to **1–2 % for N′ = 5…16** (both parities). The
    super-Schwarzian (LMRS) overshoots by 15–60 % over the same range.
  - The residual has a small systematic drift (−2 % at N′ = 7 to +2 % at N′ = 16).
  - In the j = ±1/3 sectors the agreement degrades (−2 % → +7 %, N′ = 8 → 14).
  - This confirms the D1 diagnosis ("no conformal window"). At N′ ≤ 16 we are in a double-scaled regime
    (λ ≈ 1.1–3.6), where the chord theory, as the UV completion, describes the crossover the Schwarzian misses.
- **Caveat.** Since DSSYK's small-λ limit is not the fixed-p large-N theory, the agreement must eventually fail at
  large N′ (where LMRS should take over). The growing drift may be the start of that.
- **Action item 2 (gravity computation of X and e₃).** BLY §5 is the natural route. Build the one-particle
  wormhole with the site-local crossing rules (matter–matter weight 1 for disjoint sites), solve for its zero-energy
  state, and evaluate (5.18). That is a well-posed, finite-λ computation of the zero-temperature OTOC/TOC. Its
  output can be compared directly with the measured X ≈ 0.73–0.85 and e₃ ≈ 2.3–2.4
  (`research/notes/chord_invariants_2026-09-28.md`).
- **Supersedes an old claim.** The project's older "chord diagrams with Narayana counting" gravity picture (CLAUDE.md
  overview) is not what the actual N = 2 chord theory does. The real chord rules carry arrows, friend/enemy weights,
  and matter weights q^Δ, and their planar limit is not a Haar/free-probability model.
