# Berkooz, Brukner, Narovlansky & Raz (2020) — The double scaled limit of super-symmetric SYK models

**Citation.** M. Berkooz, N. Brukner, V. Narovlansky, A. Raz, *The double scaled limit of Super-Symmetric SYK
models*, JHEP 12 (2020) 110, arXiv:2003.04405. File: `papers/berkooz_brukner_narovlansky_raz_2020.pdf` (75 pp).
("BBNR" in project documents; it is reference [1] of Boruch–Lin–Yan.)

**Note written 2026-10-05.** §§1–9 read in full. App. A read; App. B–E (inner-product derivation, eigenvector
overlaps, ground-state density resummation, conformal limit of the 2-pt function) skimmed.

## Main question
Solve the N = 2 SUSY SYK model in the double-scaled limit: density of states (fixed charge and with chemical
potential), number of ground states, 2-point functions, and the algebraic structure of the chord transfer matrix.

## Physical system
- One-flavor N = 2 SYK: N complex fermions, Q = Σ_J C_J Ψ_J with |J| = p odd, H = ½{Q, Q†} (eq. 3.4).
  - Note the ½ relative to BLY and to our code.
  - Variance ⟨C_J C_J'*⟩ = binom(N,p)^{-1} 2p J² δ (eq. 3.3). Normalized trace tr 1 = 1, so tr H = 1.
- U(1)_R generator γ = (1/2p) Σ (ψ̄ψ − ψψ̄) (eq. 3.5). The charge sector is s = (m − N/2)/p for m particles.
  Our half filling is s = 0.
- Double scaling (3.7): N, p → ∞, λ = 2p²/N fixed, q = e^{−λ}.

## Important definitions
- **Oriented chord diagrams (X–O diagrams).**
  - Each Wick pair C_J C_J* is a chord from O (a Ψ_J, i.e. a Q) to X (a Ψ̄_J, i.e. a Q†).
  - Moments are sums over k-chord diagrams (eqs. 4.1–4.3).
- **Friends and enemies (eq. 4.4).** For two index sets, the configurations ...Ψ_I...Ψ_J...Ψ̄_I...Ψ̄_J... and
  ...Ψ_I...Ψ_J...Ψ̄_J...Ψ̄_I... force I ∩ J = ∅ ("enemies"). The configuration ...Ψ_I...Ψ̄_I...Ψ_J...Ψ̄_J... allows
  overlaps ("friends").
- **Exact disentangling identity (eq. 4.7), the key finite-N statement.** Assuming no triple intersections of index
  sets, tr(Ψ_{J1}Ψ̄...Ψ_{Jk}Ψ̄...) = (−1)^{#int} 2^{−|J1∪…∪Jk|}, and |∪J| = kp − Σ_{i<j} m_ij.
  - Each friend overlap therefore contributes a factor 2^{m_ij}.
  - Averaging over Poisson(λ/2) overlaps gives friends q^{−1/2} and enemies q^{1/2} (eqs. 4.8–4.10).
- **Chord partition function (4.10):** m_k = 2^{−k} Σ_π (−1)^{#int} q^{(#e − #f)/2}.
- **Locality trick (§4.3).**
  - Two of the twelve pair configurations are non-local.
  - BBNR rewrite them via the fixed-(k, m) relations (4.11)–(4.14), where m = n_→ − n_← is a net orientation.
  - Fourier transforming m gives a local transfer matrix in a fixed charge sector, T_s = T(iλs): opening an
    O-chord costs q^{s}, opening an X-chord q^{−s} (eqs. 4.25–4.27).
  - This is the origin of BLY's m-label.
- **Inner product (§4.6, App. B).**
  - The auxiliary Fock space of X/O strings carries a *mixed q_ij-Gaussian* (Bożejko–Speicher) inner product
    (4.33–4.34), with q_ij = (q − 1)δ_ij − q (eq. 4.32; derived via eq. 8.10).
  - States with two adjacent equal chords are null.
  - The physical Hilbert space is H_phys = {|∅⟩, |n,X⟩, |n,O⟩, |n−½,X⟩, |n−½,O⟩} (eq. 4.36).
  - Its Gram matrix is given in eq. (4.37).
- **Remark after (4.34) [lit].** Several SYK operators with different λ_i give generalized statistics
  q_ij = e^{−√(λ_iλ_j)}. Tensor products give q_ij = exp(−Σ_a α_i^(a)α_j^(a)). This is the multi-species
  (Pluma–Speicher, their ref. [41]) structure.

## Main assumptions
- Double scaling. No triple intersections (corrections O(1/N)). Index overlaps are Poisson(p²/N).
- Gaussian couplings, though pairing dominance also holds more generally ([38]).
- The ground-state counting (§5.3) is a heuristic for generic couplings ("no full proof"). It may miss O(1) states
  at s = ±1/2.

## Main analytical results
1. **Transfer matrix on H_phys:** eqs. (4.38)–(4.42), with Q and Q† acting explicitly.
2. **Spectrum (§5.1–5.2).**
   - T_s is diagonalized on the exact subspaces B = Q F, B̄ = Q† F via continuous q-Hermite polynomials in q²
     (eqs. 5.4–5.10).
   - Eigenvalues: Λ_s(φ) = 2q^{−1/2}[cosh(λs − λ/2) − cos φ] (5.3). Energies: E(s,φ) = q^{−1/4}(cosh λ(s+½) − cos φ)
     (5.14).
   - Continuous density at fixed charge: eq. (5.18). Minimum energy E_min(s) ≈ (λ²/2)(|s| − ½)² (5.19).
   - With chemical potential: (5.15).
3. **Ground states (§5.3).**
   - Exact finite-N cohomology count D(r) = Σ_n (−1)^n binom(N, r + np) (5.25), nonzero only for |s| < ½.
   - Double-scaled limit: D(s) = ϑ₂(πs, e^{−π²/λ}) (5.28).
   - Small-λ form: 2e^{−π²/4λ} cos πs (5.30).
   - **For p = 3 the exact count is (2/3) cos(πs) 3^{N/2}** (text after 5.32).
4. **Super-Schwarzian limit:** the density reduces to the N = 2 super-Schwarzian (5.36–5.41).
5. **Super-Liouville:** T_s → λ²[e^ϕ − 4∂_ϕ² + (s ± ½)²] on the two bosonic subspaces (6.4–6.6), plus a fermion zero
   mode.
6. **2-point functions (§7).**
   - *Single-chord* operators A = Σ C̃_I Ψ_I (|I| = p′, charged) are computed in closed form (7.15–7.20).
   - The ground states contribute separate terms I_1 that are **not** of conformal form. So "the conformal ansatz is
     inconsistent" (7.22).
7. **Double-chord operators (§7.2).** These are O = Σ O_{IJ} Ψ_I Ψ̄_J.
   - **A Q-chord that does not cross a neutral double chord is friends with one leg and enemies with the other, so it
     is transparent. Only intersections cost: q̃^{1/2} each.** Here q̃ = e^{−2p(p′+p″)/N}.
   - So m = ⟨∅|T^{k1} q̃^{N̂} T^{k2}|∅⟩ (7.25), with N̂ the chord number operator.
8. **Quantum group (§8).**
   - The chord algebra is a contraction of sl_q(2|1) (8.10–8.15).
   - "T^{1/2}" = Q + Q† is a super-coordinate (8.16).
   - **Four-point functions are not computed.** They are listed as an open problem, with the suggestion to use
     sl_q(2|1) as in the Majorana 4-pt function of BINT [40].

## Important equations for us
(4.7), (4.9)–(4.10), (4.32)–(4.34), (4.37), (5.3), (5.14), (5.18), (5.25), (5.28), (7.25).

## Numerical methods
None beyond comparison with the known p = 3 ground-state counts ([33], [42]). Figures 8–9 plot the analytic formulas.

## Limitations
- Strictly double-scaled. Finite-p effects enter through triple intersections and through the hypergeometric (not
  Poisson) overlap statistics. These are exactly the 6–10% chord-weight corrections we measured
  (`research/notes/chord_size_dictionary_2026-10-02.md` §4).
- No 4-point functions, no matter–matter crossings, no OTOC.
- The ground-state count is heuristic at s = ±½.

## Relationship to our project
1. **The BPS dimension we observed, d = 2·3^{N′/2−1} at p = 3, half filling, is BBNR's (2/3)cos(πs)3^{N/2} at s = 0.**
   It is the cohomology count (5.25). It can now be cited rather than called "observed".
2. **Two complementary descriptions of the double-scaled theory.**
   - The fixed-charge density (5.18) is what our Krylov statement P_n = a_λ p_n^DS(0)² needs. The orthonormal
     polynomials of this measure at E = 0 should reproduce BLY's P_n. This is a check we can run.
   - BLY's Hartle–Hawking state lives in the same Hilbert space (their §2 refines BBNR §4.6–4.7).
3. **The exact trace identity (4.7) is valid at finite N whenever there are no triple intersections.**
   - At p = 3, N′ ≤ 14 triple intersections are not negligible. That is one microscopic source of the finite-p
     deviations.
   - A finite-N chord computation must use exact overlap statistics (hypergeometric, triple intersections) rather
     than q-factors.
4. **The mixed-q_ij Fock space (4.33–4.34) is the natural home for the four-arc problem (T₄, R).**
   - Our decoder has four chord species: O, X for Q and O′, X′ for Q′ = UQU, or equivalently for an OU-correlated
     copy.
   - With independent copies (x = 0) this is the Pluma–Speicher multi-flavor construction. Cross-flavor crossings
     carry q, because index sets of different copies overlap with the same Poisson law.
   - With correlated couplings, the pairing between flavors carries x, while crossings are unchanged.
   - *Corrected 2026-10-05.* An earlier version of this note called this "a correlated mixed q-Gaussian, beyond
     Bożejko–Speicher and Liu–Shen". That is wrong, at least in the Majorana case. By the multilinearity argument of
     Pluma–Speicher Cor. 3.2, correlated copies are the ordinary q-Gaussian functor evaluated on non-orthogonal
     vectors ⟨h, h_x⟩ = x. Conditioning acts as x^{N̂} (Bożejko–Kümmerer–Speicher Thm. 2.11). The N = 2 oriented-chord
     version is expected but not written down anywhere we have seen.
   - The genuinely new object is the **four-arc** function. Its adjacent-arc chords weigh x and its opposite-arc
     chords weigh 1 (a ℤ₂ rule), unlike BINT's random-matter product rule (`docs/derivations.md` D5).
5. **The double-chord rule (§7.2) is the chord rule for our probe.**
   - The single-site bilinear has p′ = p″ = 1, so q̃^{1/2} = e^{−2p/N} = x_DS per crossing.
   - The $\mathbb Z_2$ property of U = (−1)^{n_i} (crossing twice cancels) is *not* in BBNR's random-operator rule,
     where two crossings would give q̃. Our probe is a fixed operator, so the rule must be modified.
6. **Ground states spoil the conformal ansatz (§7.1.1).** This is consistent with our decoder being a pure
   ground-state (BPS-projected) quantity. Any Schwarzian/JT four-point function used for T₄ must be the
   zero-energy (LMRS-type) one, not the conformal ladder.

## Note added 2026-10-14 (D7)

- **Sign in (5.7)/(5.10).** With T b_n as in (5.5), the coefficients α_n = q^{n/2}H_n(cos φ|q²)/(q²;q²)_n give the
  eigenvalue q^{−s} + q^{s−1} + 2q^{−1/2}cos φ, i.e. +cos φ. The eigenvector for Λ_s(φ) of (5.3) therefore has
  H_n(−cos φ|q²). This was checked directly against (5.5) (`tests/test_chord_dressing.py::test_bbnr_eigenvector_sign`).
  It is harmless for BBNR's density of states (the measure is symmetric under φ → π − φ). It matters for any kernel with
  one zero-energy end: combining (5.10) with (5.3) reflects the BPS→continuum weight across the band.
- **Used here.** §7.1 (charge shift s → s ± p′/2p, the q̃^N̂ junction factor) is the matter rule behind our closed-form
  BPS→continuum two-point function (D7 §1). Our transfer-matrix engine (`src/chord_transfer.py`) implements the chord
  rules (4.16)–(4.17) and the inner product (4.37), and reproduces (4.40)–(4.42) and (5.5). (7.17) contains the
  ground-state/continuum moments integrated over the charge; D7 gives the fixed-charge spectral measure in closed form.
