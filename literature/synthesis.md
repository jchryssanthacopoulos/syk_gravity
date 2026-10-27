# Literature synthesis

*Last revised 2026-10-27 (sixth pass: added the LMRS companion 2207.00407, Penington–Witten 2024 and Lee–Lee 2019; §1 block 2026-10-27, §3 items 14–15, §6 updated, with a correction to the earlier Lee–Lee entry). Fifth pass 2026-10-26 ( added Kanazawa–Wettig 2017, the per-degree Hodge/RMT classification; §3 item 13).
Fourth pass 2026-10-05 ( added BBNR 2020, BINT 2018, Pluma–Speicher 2019, Bożejko–Kümmerer–Speicher
1997, G. Lin 2026; new §8 on tools for the T₄/R program). Third pass 2026-09-30 added Belin–Fu–La Rocca 2026, Lin 2022,
Liu–Shen 2026 and §7 on where the project sits. Second pass 2026-09-28 added 9 papers: PSSY, Turiaci–Witten, Boruch–Iliesiu–Yan, HITZ,
Heydeman–Turiaci–Zhao, FGMS, Y. Chen, Collins–Matsumoto–Novak, Berkooz–Mamroud. Per-paper notes are in
`literature/notes/`; BibTeX in `literature/bibliography.bib`.*

## 1. The one-paragraph picture

Fortuitous BPS states (Chang–Lin 2024; Chang–Chen–Sia–Yang 2024; the N=2 SYK model of FGMS 2016) are the
candidate black-hole microstates. They are expected to be *strongly chaotic*. Because their energies are exactly
degenerate, chaos must be diagnosed through **LMRS operators** Ô = P_BPS O P_BPS (LMRS 2022; Chen–Lin–Shenker
2024), Berry curvature (Chen–Colin-Ellerin–Mamroud–Papadodimas 2026), or random-matrix models of the spectrum
(Turiaci–Witten 2023; Johnson 2026). On the gravity side, N=2 JT supergravity is dual to independent AZ ensembles
per R-charge (Turiaci–Witten). BPS states at fixed R-charge number e^{S0}cos(πQ) (HITZ App. A). States prepared by
matter insertions behave as **Haar-random vectors inside the BPS space**, with Marchenko–Pastur Gram statistics and
*no* higher-genus corrections (Boruch–Iliesiu–Yan, the BPS analogue of PSSY). Finite-Δ crossing weights are
supplied by chord diagrams (Berkooz–Mamroud). Our operator O_T = P_B Π_T P_B (Chryssanthacopoulos–Vegh 2026) is an
LMRS operator whose UV operator Π_T = 1 − n_{N′} is a *light* fermion bilinear (Δ = 1/q̂ = 1/3 from FGMS). It is also
canonically tied to uplifting BPS states from N to N+1. Its spectrum has so far been compared against
**Wachter/MANOVA**, the law of two free projections. That law is realized exactly by Haar/Jacobi ensembles
(Collins 2004), has a 1/D² expansion given by signed monotone-walk counts (Collins–Matsumoto–Novak 2021), and is
studied for structured projections (Haikin–Zamir–Gavish 2017; Kunisky 2023) and via bordered-Toeplitz data
(Dubbs–Edelman 2015) and topological recursion (Bouchard 2024).

*Added 2026-09-30.*
- Belin–Fu–La Rocca 2026 make the BPS-chaos hypothesis "the BPS subspace is Haar-oriented" precise. For projected
  operators it gives GUE/semicircle at small relative dimension, the Wachter law for a spin-½ probe, and a generic
  count of wall atoms. It is exactly our Haar null model.
- Lin 2022 shows that in double-scaled SYK the **chord number equals operator size** (divided by the size of one H
  term). Two-sided light-probe correlators measure it.
- Liu–Shen 2026 prove that independent or overlapping SYK Hamiltonians converge to mixed q-Gaussians, with
  inter-family crossing weight q_ij = (−1)^{r_i r_j}e^{−2λ_ij}. Freeness needs λ_ij → ∞.
- Together these place our exact finite-N results at the junction of three literatures: fortuity/BPS chaos, DSSYK
  chords/operator size, and free probability (§7).

*Added 2026-10-05.*
- BBNR 2020 solve N=2 SUSY SYK in double scaling: oriented chords with friend/enemy weights, fixed-charge densities,
  and the exact finite-N ground-state count. That count explains our d = 2·3^{N′/2−1}.
- BINT 2018 give the matter-chord rule q̃ = e^{−2pp̃/N} and the crossed 4-pt function / R-matrix of Majorana DSSYK.
- Pluma–Speicher 2019 and Bożejko–Kümmerer–Speicher 1997 make "correlated couplings ⇒ x per chord ⇒ x^{N̂} on chord
  space" rigorous (Majorana, double scaled). They also show our four-arc OU function is an **out-of-time-order
  correlator of the q-OU process**.
- G. Lin 2026 gives the first explicit zero-energy N=2 super-JT OTOC (equal Δ). His Berry-curvature second moment is
  TOC − OTOC.
- See §8 for what this means for the T₄/R program.

*Added 2026-10-26.*
- Kanazawa–Wettig 2017 (KW) classify N=2 SYK level statistics per fermion number. They do this per *Hodge piece*:
  V⁺ (coexact), V⁻ (exact) and V^z (BPS), with tables of all three dimensions for N_c = 3–17 (App. B).
- The cohomology premise is CV's fortuity window.
- They explain the exceptional window-edge zero modes at N_c ≡ 1 (mod 4) by Kramers degeneracy, but not the excess at
  N_c = 9.
- Every V^± is GUE except one supermultiplet: the one that straddles the window, from coexact at ℓ = −q̂ to exact at
  ℓ = +q̂. It is GOE or GSE.
- For the uplift anywhere in the window (D10), that multiplet is where off-centre dressing resonates.

*Added 2026-10-27.*
- **Gap per R-charge (LMRS companion eq. 33; Penington–Witten 4.133–4.135, 4.156).** Outside the BPS window the
  continuum starts at E_gap(j) = (|j| − ½)²/(8C), with C = α_S J N. LMRS check it in N=2 SYK at N = 16 (their Fig. 9).
  PW derive the same threshold structure in canonical super-JT; it closes at |J| = ½.
- **[derived from (33); interpretation of our D10 regimes]** A dressing K = Q₀⁺E passes through the supermultiplet
  coexact(ℓ′−1) ↔ exact(ℓ′+2q−1). Its two charges are j and j + 1, with j + 1 = s_in (Q carries R-charge 1 when
  fermions carry 1/q̂).
  - For j ∈ (−1, 0), (33) gives the same gap at both charges, (|j| − ½)² = (|j+1| − ½)² = (s_in − ½)².
  - So the lowest energy of the intermediate multiplet is (s_in − ½)²/(8C).
  - **gapped** (s_in > ½): finite threshold, so Q₀⁺ is bounded.
  - **resonant** (s_in = ½): the threshold closes, which matches the heavy K^{−1/2} tail found in D11 layer 1.
  - **obstructed** (s_in < ½): the threshold is again finite. The obstruction there comes from BPS states at the target
    charge, not from gaplessness.
  - Status: the algebra is a derivation from (33). That this threshold sets the dressing norms is an interpretation.
    A cheap test: at fixed N′, the smallest nonzero singular value² of Q₀ on that multiplet in the three gapped h
    sectors ℓ′ = −1, 0, 1 (q = 3: s_in = 2/3, 5/6, 1) should approach the ratio 1 : 4 : 9. The curved CCSY edge
    (D7 [conj]) may shift this at small N.
- **Penington–Witten.**
  - Exactly one simultaneous zero-energy state per SU(1,1|1) matter multiplet with |2u − J_R| < ½, BPS *or long*.
  - The projected algebra is Type II₁. Neutral projected operators multiply in closed form (5.20–5.21).
  - The trace weights R-charge sectors by cos(πJ_R) (5.62), i.e. the BPS counting.
  - The projected two-point function 2^{−4Δ}Γ(2Δ+1) (5.32) equals LMRS-short (10): an independent check of the
    normalization.
  - Π_T is a composite (1 − n_{N′}), not a primary, so PW do not directly give the decoder's moments.
- **Lee–Lee.** The N=2 global multiplets, i.e. su(1,1|1) lowest-weight representations: L[j₀]_h with h > |j₀|/2, and
  A[j₀]_{|j₀|/2}.
  - Recombination L → A ⊕ A′ (4.4); characters (6.10–6.12).
  - Dictionary to PW, derived here: λ = h₀, u = j₀/2.
  - They do **not** contain OPE, tensor-product or exchange data. The super-6j input of G. Lin 2026's OTOC route (his
    eq. 6.8) is still missing from the literature we hold.

## 2. Literature map

```
 FGMS 2016 (N=2 SYK, index) ──► HTZ 2022 (phases; two-flavor ψψχ) ──┐
          │                                                          │
 Chang–Yin 2013 ─► Chang–Lin 2022 ─► Chang–Lin 2024 (fortuity) ──► CCSY 2024 (fortuity in SYK; monotone tower)
                                                                     │        Chen 2025 (single matrix; solvable)
                                                                     ▼            │
 ─── gravity / BPS chaos ───────────────────────────── Chryssanthacopoulos–Vegh 2026 (decoder D_T, O_T; ⟨r⟩=GUE)
 PSSY 2019 (EOW branes, Gram/MP) ─► BIY 2023 (BPS states from matter; MP; no genus; Haar in BPS)         │
 HITZ 2020 (N=4/N=2 super-Schwarzian; e^{S0}cos πQ; gap)                                               │
 Turiaci–Witten 2023 (N=2 JT ↔ AZ ensembles per R-charge) ─► Johnson 2026 (Bessel/Airy BPS model)      ▼
 LMRS 2022 (Ô=POP; Δ-controlled crossings) ─► Chen–Lin–Shenker 2024 ─► Berry 2026          current draft
 Berkooz–Mamroud 2024 (chords: crossing weights q^{ΔΔ'}; SUSY DSSYK)   Miyahara–Shibuya 2026   (gravitational_dual_
   └─► Lin 2022 (chord Hilbert space; chord # = size/q) ─► BLY 2023 (N=2 super-chords; P_n)      fragility)
 Belin–Fu–La Rocca 2026 (Haar-oriented BPS subspace ⇒ GUE/Wachter; generic atoms) ── BPS-chaos hypothesis = our null
 ─── random matrices / free probability ────────────────────────────────────────────────────────▲
 Collins 2004 (exact Jacobi law) ─ Collins–Matsumoto–Novak 2021 (Weingarten 1/N = monotone walks) │
 Kunisky 2023 (MANOVA theorems, error identity) ─ Haikin–Zamir–Gavish 2017 ─ Dubbs–Edelman 2015    │
 Demni–Hamdi–Hmidi 2012 (free Jacobi/liberation) ─ Liu–Chen–Balents 2017 ─ Bouchard 2024 (TR) ────┘
 Liu–Shen 2026 (SYK families → mixed q-Gaussians, q_ij = ±e^{−2λ_ij}; ε-freeness; Poisson overlap lemma)
 ─── added 2026-10-05 ──────────────────────────────────────────────────────────────────────────────
 BINT 2018 (matter chords q̃; crossed 4-pt; R-matrix = U_q(su(1,1)) 6j) ─► Berkooz–Mamroud 2024 (review)
 BBNR 2020 (N=2 DSSYK: oriented chords, friend/enemy, ground-state count, double-chord rule) ─► BLY 2023
 BKS 1997 (q-Gaussian functor Γ_q(T); q-Mehler; q-OU is Markov) ─► Pluma–Speicher 2019 (SYK copies & Brownian
   couplings → q-Gaussian/q-Brownian) ─► Liu–Shen 2026
 LMRS 2022 ─► G. Lin 2026 (particle on SU(1,1|1): exact zero-energy N=2 OTOC; Berry 2nd moment = TOC − OTOC) ◄─ CCMP 2026
 ─── added 2026-10-26 ──────────────────────────────────────────────────────────────────────────────
 FGMS 2016 ─► Kanazawa–Wettig 2017 (AZ classes per Hodge piece V^±_f; exceptional zero modes = Kramers; App. B dims) ─► CV 2026 [25]
 ─── added 2026-10-27 ──────────────────────────────────────────────────────────────────────────────
 HTZ 2022 (gap) ─► LMRS-short 2022 (E_gap(j) = (|j|−½)²/8C, N=16 check; Γ(1+2Δ)/4^{2Δ}; fuzzball depth log S_e) ─► window edge = resonance (D10/D11)
 LMRS 2022 ─► Penington–Witten 2024 (one E_L=E_R=0 state per SU(1,1|1) multiplet; Type II₁ projected algebra; cos πJ trace) ─► G. Lin 2026 [49]
 Lee–Lee 2019 (su(1,1|1) multiplets, recombination, characters; no exchange data) ─► G. Lin 2026 §6 [61]
```

| Paper | Provides | Used by the project for | Connects to |
|---|---|---|---|
| Chryssanthacopoulos–Vegh 2026 | decoder D_T, K_T, O_T; fidelity/dressing cost; LES; Table 1; ⟨r⟩ = 0.599 | everything | CCSY (model), LMRS (operator type) |
| FGMS 2016 | N=2 SYK; Δ_ψ = 1/(2q̂); twisted index; exact BPS degeneracies (eq. 5.7) | model, BPS dims, Δ of Π_T | CCSY, HTZ, Turiaci–Witten |
| HTZ 2022 | charge-dependent Δ; phase transition; two-flavor ψψχ model | two-flavor control | CCSY, CV |
| CCSY 2024 | fortuity in SYK; R-charge concentration; symmetrized two-flavor monotone tower; LMRS test | model, BPS counts, control | CV, Chang–Lin |
| Chen 2025 | solvable single-matrix model with fortuity, no chaos | CV §6.1 control | CCSY |
| Chang–Lin 2024 / 2022, Chang–Yin 2013 | fortuity concept; 1/16-BPS cohomology history | background | — |
| LMRS 2022 | super-JT zero-energy correlators of Ô; Δ-controlled crossings; N_BPS = e^{S0}L̂ | the gravity calculation to use (Q1, Q4) | Chen–Lin–Shenker, Berry, BIY |
| Chen–Lin–Shenker 2024 | LMRS criterion; strong vs weak chaos; width ~ S0^{−Δ} | interpretation of ⟨r⟩ | LMRS |
| Berry 2026 | chord technology at Δ→∞; free commutator; SYK BPS numerics; UV imprint | draft's chord section | LMRS, Berkooz–Mamroud |
| Miyahara–Shibuya 2026 | BPS-subspace chaos↔integrability crossover with fortuitous states | test of CV conjecture (Q6) | CCSY |
| PSSY 2019 | EOW-brane states; Gram-matrix MP from planar replica wormholes; Gaussian-state dual | cited as template by draft | BIY |
| BIY 2023 | BPS states from matter insertions; MP Gram; genus/trumpets vanish; Haar-random states in H_BPS | correct BPS analogue of draft's model (Q4) | PSSY, LMRS, HITZ |
| HITZ 2020 | super-Schwarzian spectra; ρ_ext = e^{S0}cos(πQ); gap | BPS normalization per sector | Turiaci–Witten, LMRS |
| Turiaci–Witten 2023 | N=2 JT ↔ AZ(1+2ν,2) ensembles per R-charge; fixed BPS counts; SYK check | energy-side RMT | Johnson, CCSY |
| Johnson 2026 | energy-side Wishart/Bessel model; repro scripts | contrast only | Turiaci–Witten |
| Berkooz–Mamroud 2024 | chord rules incl. matter (q̃ = q^Δ); SUSY DSSYK; 1/N^p fluctuations | finite-Δ crossing weights (Q1, Q4) | LMRS, Berry |
| Kunisky 2023 | freeness criterion for ABA; exact error identity | draft's "error formula" | Collins, HZG |
| Collins 2004 | exact finite-n Jacobi law; edge universality; large deviations | Haar null model | Kunisky |
| Collins–Matsumoto–Novak 2021 | Weingarten calculus; 1/N expansion = signed monotone walks | Haar null model genus expansion | Collins |
| Dubbs–Edelman 2015 | Jacobi parameters of Wachter; bordered-Toeplitz moment problem | draft's ρ_dec | Kunisky |
| Haikin–Zamir–Gavish 2017 | structured frames have MANOVA spectra | caution | Kunisky |
| Demni–Hamdi–Hmidi 2012 | free Jacobi (liberation) process | model for partial freeness (Q2) | Collins |
| Liu–Chen–Balents 2017 | Wachter law in SYK₂ entanglement | draft's EE "bridge" | — |
| Bouchard 2024 | topological recursion conventions | draft's W₁ | Turiaci–Witten (loop eqs.) |
| Belin–Fu–La Rocca 2026 | projected spectrum under Haar orientation: exact B-spline density (2.33); GUE/semicircle for α ≪ 1 (2.65–2.66); spin-½ ⇒ Wachter-type law (3.3) and eigenvalue recovery (3.8); generic atom count (Thms. A.1–A.2); free compression (4.69) | the Haar null model stated as the BPS-chaos hypothesis; generic atom count vs SYK atoms (Q3); open question "how Haar is the physical rotation?" | LMRS, Chen–Lin–Shenker, Kunisky, Collins |
| Lin 2022 | DSSYK chord Hilbert space = bulk Hilbert space; chord number = Krylov (46–48); two-sided light-probe correlators measure length (53); **chord number = size/q** (57) with O(λ) smearing (59) | origin of size ↔ wormhole-length; our variance formula = exact analogue of eq. 53; D3 (size per chord unit = p−1) | BLY (N=2), Berkooz–Mamroud |
| Liu–Shen 2026 | SYK with different lengths / overlapping supports → mixed q-Gaussians, q_ij = (−1)^{r_ir_j}e^{−2λ_ij} (Thms. 1.1–1.2); Poisson overlap lemma (4.5) and its failure for r/n ↛ 0 (Ex. 4.6); asymptotic ε-freeness | rigorous inter-family crossing weights; two independent SYK are not free at finite λ; our finite-N factorization caveat and x_s = 0 exactness | Berkooz–Mamroud, Speicher |
| BBNR 2020 | N=2 DSSYK: oriented X–O chords, friend/enemy weights q^{∓1/2} (4.4–4.10); exact trace identity without triple intersections (4.7); fixed-charge transfer matrix (4.25–4.27); mixed-q_ij inner product (4.32–4.37); fixed-charge density (5.18); **finite-N ground-state count (5.25), = (2/3)cos(πs)3^{N/2} at p = 3**; double-chord (neutral bilinear) rule q̃^{N̂} (7.25); sl_q(2\|1) (§8); 4-pt functions listed as open | citable origin of d; Krylov check via (5.18); home of the four-chord-species (Q, Q†, Q′, Q′†) problem; probe rule | BLY, BINT, Berkooz–Mamroud |
| BINT 2018 | matter chords q̃_A = e^{−2pp_A/N} (2.9–2.11); 2-pt vertex (5.4) = q-Mehler kernel; bi-local conserves energy (3.20); uncrossed (4.9) and crossed (4.13) 4-pt; R-matrix (5.6), integral form (6.3), Schwarzian 6j limit; U_{q^{1/2}}(su(1,1)) (§7); finite-p structure (§8) | the crossed-channel template; shows the random-matter product rule differs from our ℤ₂ four-arc rule (D5) | BBNR, Berkooz–Mamroud, Lin 2022 |
| Pluma–Speicher 2019 | iid SYK copies → q-Gaussian family (Thm. 3.1); Brownian couplings → q-Brownian motion by multilinearity (Cor. 3.2); exact word sign (Lemma 3.3); higher cumulants at scale binom(n,p)^{−1} (Thm. 4.1); no operator-valued continued fraction (§6) | OU-correlated copies = q-Gaussians with ⟨h,h_x⟩ = x (corrects our BBNR note); fluctuation scale of per-realization quantities | Liu–Shen, BKS |
| BKS 1997 | q-Hermite & Rogers q-Mehler kernel (Thm. 1.10); second quantization Γ_q(T) for contractions, CP and trace-preserving (Thm. 2.11); q-OU is the unique stationary Markov q-Gaussian (Cor. 3.10, Prop. 3.11); classical versions reproduce only time-ordered correlators (§4) | rigorous "x per chord = x^{N̂}"; T₄^OU is out of OU-time order, so it is not fixed by the Markov kernel (D5) | Pluma–Speicher, BINT |
| G. Lin 2026 | super-JT as a constrained particle on SU(1,1\|1), PSU(1,1\|2); physical supercharges; Wilson lines; zero-energy 3-pt; **exact zero-energy N=2 OTOC for equal Δ (5.26–5.38)**; N=4 large-Δ (5.41–5.42); Berry 2nd moment = TOC − OTOC (6.15–6.19) | triple-scaled target for the crossed channel of our light probe (Δ = 1/3); the isotropic OU expansion sees only the quantum metric, not curvature (D5) | LMRS, CCMP, BLY |
| Kanazawa–Wettig 2017 | AZ classes of N=0,1,2 SYK. For N=2: V_f = V⁺_f (coexact) ⊕ V⁻_f (exact) ⊕ V^z_f (BPS) (6.8–6.10). Pairings N⁺_f = N⁻_{f+3} (6.12), N^±_f = N^∓_{N−f} (6.11). Zero modes in abs(f − N_c/2) ≤ 3/2 (premise 6.15). Exact sequences and closed forms for N^± (6.16–6.21; Table 8). Exceptional zero modes via DIII Kramers (§6.2; N_c = 9 excess unexplained, fn. 12). GOE/GSE only for the multiplet V⁺_{(N_c−q̂)/2} ↔ V⁻_{(N_c+q̂)/2} (6.13, 6.24). App. B tables. | Hodge dimensions at every degree, used by the uplift anywhere in the window (D10); the exceptional-state counts (verified against our code); the window-edge multiplet that controls resonant dressing | CV ([25]), FGMS, Miyahara–Shibuya, BBNR |
| LMRS-short 2022 (2207.00407) | Zero-energy 2-pt e^{S0}Γ(1+2Δ)/4^{2Δ} (10); ⟨ℓ⟩ = 2(γ_E + log 4) (11); c_Δ² ∝ S_e^{−2Δ} (24), depth d̃ = log S_e (26); disk = cylinder in the probe limit (27–29); large-Δ POP semicircle; gap per R-charge (33) checked at N=16 (Fig. 9); Table 1 SYK vs Schwarzian | Gravitational origin of the window-edge resonance (D10/D11); normalization of projected correlators | LMRS 2022, HTZ, PW |
| Penington–Witten 2024 | Canonical N=2 super-JT plus matter: one simultaneous ground state per SU(1,1\|1) multiplet with abs(2u−J_R)<½, even non-BPS (§4.4); thresholds (4.133–4.135); Type II₁ projected algebra with closed product (5.20–5.21); cos πJ_R trace = BPS counts (5.62); App. E highest-weight gauge ↔ chord space (bosonic) | Exact zero-energy count per multiplet (uplift dictionary, C1); algebraic status of the decoder; route to a gravity derivation of chord rules | LMRS, G. Lin 2026, BBNR, HITZ |
| Lee–Lee 2019 | 2D global superconformal multiplets for all N; N=2 Table 2 (h ≥ abs(j₀)/2); recombination (4.4); characters (6.10–6.12) | Multiplet labels and channel sums for an N=2 OTOC (C3/A3); recombination as an algebraic analogue of fragile pairs (speculative) | G. Lin 2026, PW |

## 3. Established results (literature) that the project relies on

1. **Model and BPS counts.** One-flavor N=2 SYK (FGMS eq. 5.1). Exact BPS degeneracies for q̂=3 (FGMS eq. 5.7;
   Berry Table 2; CV Table 1): N even: 2·3^{N/2−1} at Q_R = 0, 3^{N/2−1} at ±1/3; N odd: 3^{(N−1)/2} at ±1/6
   (+1 or 3 at ±1/2 if N ≡ 1 mod 4). Entropy density ½ln3 from the twisted index (FGMS eqs. 5.5–5.6).
   Concentration in the fortuity window (CCSY eqs. 1.8–1.9).
2. **Decoder.** D_T, K_T, O_T and their relations (CV eqs. 42–63, 101–102); GUE adjacent-gap statistics of O↓ at
   N′=12, P=6 (CV eq. 105).
3. **Conformal data.** Δ_ψ = 1/(2q̂) at zero background charge (FGMS), hence Δ = 1/3 for the bilinear n_i − ⟨n_i⟩
   at q̂=3. Away from zero charge Δ depends on the charge (HTZ eq. 2.27).
4. **N=2 JT / super-Schwarzian.** BPS density per R-charge e^{S0}cos(πQ), |Q| < 1/2, plus a gap
   (HITZ eq. A.18; LMRS eq. 38). Dual ensemble: independent AZ(1+2ν,2) blocks per R-charge with non-fluctuating
   BPS counts (Turiaci–Witten §2.3).
5. **BPS-projected operators.** Zero-energy correlators are nontrivial functions of Δ and R-charge (LMRS
   eqs. 65–66). Disk crossings are suppressed only at large Δ (LMRS eq. 159; Berry p. 60). Chord crossings weigh
   q^{Δ}-type factors, → 1 for fixed-size operators at large N (Berkooz–Mamroud eq. 2.28).
6. **BPS state preparation.** Matter-prepared BPS states have MP Gram statistics with dimension e^{S_j} =
   e^{S0}cos(πj). Higher-genus and trumpet contributions vanish at β→∞. The statistics are reproduced exactly by
   Gaussian (Haar) random vectors in H_BPS (BIY eqs. 3.43–3.51, 3.68–3.72). Non-BPS analogue: PSSY.
7. **Wachter / Haar null model.** Limit of PQP for free projections (Collins §3.2; Kunisky Thm. 1.4); exact Jacobi
   law at finite size (Collins Thm. 2.2); 1/D² corrections = signed monotone-walk counts (Collins–Matsumoto–Novak
   Thm. 4.5); freeness criterion and exact error identity (Kunisky Thm. 1.5, Lemmas 4.3, 4.6).
   - Same model from the BPS-chaos side: Belin–Fu–La Rocca eqs. 3.3, 3.8, 4.69.
   - Under generic (Haar) orientation, an eigenvalue of multiplicity γ appears in the compression exactly
     max(0, K − N + γ) times, almost surely (Belin–Fu–La Rocca Thms. A.1–A.2).
8. **Chord number = operator size (DSSYK).**
   - n̄ = size/q, where q is the size of one H term (Lin eq. 57; derivation eq. 58; smearing by Poisson overlaps,
     eq. 59).
   - A two-sided correlator of a light s-fermion probe, averaged over index sets, gives e^{−Δℓ} (Lin eq. 53); s = 1 is
     the size operator (Lin eq. 55).
   - For N=2: n_O, n_X ↔ size_O/p, size_X/p (BLY eqs. 5.20–5.23).
9. **Inter-family crossings.**
   - Independent SYK Hamiltonians with λ_ij = lim r_ir_j|A_i∩A_j|/n² converge to a mixed q-Gaussian system with
     q_ij = (−1)^{r_ir_j}e^{−2λ_ij} (Liu–Shen Thms. 1.1–1.2).
   - Two independent SYK models on the same fermions cross each other with the same q as themselves, so they are
     not free at finite λ.
   - Factorization of sign averages over several chords is a Poisson limit (Liu–Shen Lemma 4.5). It fails when a set
     size is a finite fraction of n (Ex. 4.6).
10. **N=2 DSSYK ground states (BBNR).**
    - The exact finite-N cohomology count is D(r) = Σ_n (−1)^n binom(N, r + np) (BBNR 5.25), = (2/3)cos(πs)3^{N/2} for
      p = 3. Its double-scaled limit is ϑ₂ (5.28).
    - The fixed-charge continuous density is (5.18).
    - Ground states spoil the conformal ansatz for 2-pt functions (§7.1.1).
11. **Correlated couplings (Majorana, double scaled).**
    - Copies with coupling correlation x converge to q-Gaussians with ⟨h, h_x⟩ = x (Pluma–Speicher Cor. 3.2 argument).
    - Conditioning acts as x^{N̂} on chord space (BKS Thm. 2.11), so two-arc functions are q-Mehler integrals (BKS
      Thm. 1.10). These are the same kernels as BINT's matter vertex.
    - The q-OU process is Markov, but four-arc alternating functions are out of time order (BKS §4).
12. **Zero-energy N=2 crossed 4-pt function in super-JT** (equal Δ): G. Lin 2026, eqs. 5.26–5.38. The large-Δ ratio
    OTOC/TOC ~ e^{−3.52Δ} is LMRS eq. 159. For chiral operators the long-time TOC = OTOC up to sign (LMRS eq. 86 and
    Fig. 9).

13. **Hodge pieces per degree and their statistics (Kanazawa–Wettig).**
    - Exact, coexact and BPS dimensions at every fermion number, given the window premise (KW 6.16–6.21, App. B).
      Our rank computations reproduce the App. B rows for N_c = 5, 9, 12, 13.
    - Labels: KW's V⁺ is coexact and V⁻ exact. CV App. A uses the opposite labels (literature/notes/kanazawa_wettig_2017.md).
    - Level statistics on each V^±_f are GUE except on the window-straddling multiplet, which is GOE (N_c ≡ 3 mod 4) or
      GSE (N_c ≡ 1 mod 4) at q̂ = 3 (KW 6.13); general odd q̂ in 6.24.
    - That multiplet is Turiaci–Witten's "CT-invariant multiplet (−q̂/2, q̂/2)", β = 1 at N = 11 (TW Fig. 4, pp. 18–19).
      This resolves the β = 1 remark in `literature/notes/turiaci_witten_2023.md`.
    - Exceptional zero modes at the window edges for N_c ≡ 1 (mod 4) at q̂ = 3: 1, 3, 1, 1 at N_c = 5, 9, 13, 17. The
      N_c = 9 excess is unexplained.
    - For q̂ = 5, KW report edge modes only at N_c = 7, 11 up to N_c = 17. Parity forces them again at N_c = 19
      (project, D10 §3).

14. **Gap per R-charge (LMRS-short eq. 33; HTZ; PW 4.133–4.135).** E_gap(j) = (|j| − ½)²/(8C) for |j| > ½, with C = α_S J N
    and α_S = 0.00842 for q̂ = 3. Checked at N = 16 for one realization (LMRS-short Fig. 9). It vanishes at the window
    edge |j| = ½.
15. **Zero-energy states and the projected algebra (Penington–Witten).**
    - One simultaneous ground state per SU(1,1|1) matter multiplet with |2u − J_R| < ½ (§4.4).
    - Projected neutral operators close under multiplication (5.20–5.21), and the matter vacuum is the unique trace.
      The algebra is Type II₁, assuming it is a factor.
    - With charged operators, Tr Π_{J_R} ∝ cos(πJ_R) (5.62).

## 4. Tensions between the literature and the current draft

| Draft claim | What the literature says |
|---|---|
| e^{S0} = D (ambient Fock-sector dimension) | e^{S0}cos(πQ) = BPS degeneracy, i.e. ∝ d (LMRS eq. 82; HITZ eq. A.18; BIY eq. 3.15; Berry eq. 6.4) |
| BPS projection trivializes the super-Schwarzian; only e^{S0χ} survives | Zero-energy correlators depend on Δ, j (LMRS eqs. 65–66, 85; BIY eqs. 3.19, 3.43) |
| Gravity gives a handle tower ∝ D^{−2g} matching the Haar/Weingarten expansion | In the BPS sector trumpet/higher-genus contributions **vanish** (BIY §3.4). The Haar 1/D² terms are signed monotone-walk counts (CMN Thm. 4.5), not a sum over surfaces with weight D^χ |
| Crossing chord = handle ∝ 1/D² | Disk crossings are O(1) in e^{S0} and weighted by Δ (LMRS eqs. 159, 169; Berry p. 60; Berkooz–Mamroud eq. 2.28) |
| "Role of Δ played by D→∞" | No such statement anywhere. Freeness in LMRS/Berry comes from Δ→∞; Π_T is light (Δ = 1/3) |
| Slot projector Π_T as a second EOW-brane species (PSSY-style) | PSSY and BIY use one family of prepared states (EOW flavours or matter insertions). No construction of occupation-pattern branes exists |
| PSSY's MP is the b→0 limit of Wachter | Analogy at best, needs rescaling; BIY is the relevant BPS computation |
| Soft tail = image of near-BPS continuum (Johnson/HITZ density) | Tail eigenvectors are near-harmonic (audit), not near-BPS energy eigenstates |
| arXiv:2412.06902 = "Chang & Lin"; N2matrix = "Turiaci, Witten et al." | Chang–Chen–Sia–Yang; Turiaci & Witten (two authors) |
| ½ln3 via "Lagrange–Bürmann proof" | FGMS twisted index, eqs. 5.5–5.7 |
| LCB's (κ,λ) ↔ our (a,b) | LCB's λ (compressed dimension) ↔ our a; κ ↔ b |

## 5. Three different notions of "random/free" that must be kept apart

1. **Haar freeness of the BPS subspace inside the Fock space** (the draft's model; Collins, Kunisky, CMN):
   P_B random relative to Π_T ⇒ Wachter(a,b), variance ab(1−b) ∝ a → 0 exponentially in N.
2. **Haar randomness of states inside the BPS space** (BIY §3.5; Turiaci–Witten's random BPS wavefunctions;
   ETH-type): matrix elements ⟨ψ_i|O|ψ_j⟩ of a simple operator between BPS states are random with variance set
   by the zero-energy two-point function (LMRS). This is what gravity actually predicts. Its N-dependence is a
   power law, ~ S0^{−2Δ} (Chen–Lin–Shenker; LMRS eqs. 83–85).
3. **LMRS/chord freeness at large Δ** (LMRS §5.3; Berry §6; Berkooz–Mamroud): non-crossing dominance, semicircular
   Ô. Absent for light Π_T, where crossed and uncrossed diagrams have equal weight.

All three give a narrow spectrum around b at large N, but with different widths and shapes. Distinguishing (1)
from (2) through the N-scaling of Var(λ) is research question Q1. Within (2), the Gaussian-vs-free question is
controlled by Δ (3).

**Refinements (2026-09-30).**
- Belin–Fu–La Rocca argue "chaos ⇒ (1)" and derive its consequences. Our exact results quantify how (1) fails for
  SYK, which is only U(N)-invariant.
- For a U(N)-invariant projector ensemble, second-moment quantities depend only on the size spectrum w_k. That covers
  the decoder variance and the fourth moment of two independent BPS spaces.
- (1) holds at this level iff the U(N)-twirl is flat on all nonzero sizes, i.e. **freeness = flat transmission =
  vanishing 2-design defect** (`research/pdfs/progress_2026-09-29.pdf` §5).
- Notion (3) has a Hamiltonian-level rigorous counterpart: mixed q-Gaussians with inter-family crossing q_ij
  (Liu–Shen). Free ⇔ q_ij = 0 ⇔ λ_ij → ∞, i.e. heavy chords, matching LMRS/Berry's large-Δ freeness.

## 6. Gaps in the literature collection

Still missing and worth adding:
- ~~LMRS companion, arXiv:2207.00407~~ (acquired 2026-10-27; `literature/notes/lin_maldacena_rozenberg_shan_2022_short.md`).
- ~~Boruch, Lin, Yan~~ (acquired 2026-09-28).
- ~~Berkooz, Isachenkov, Narovlansky, Torrents, arXiv:1811.02584~~ (acquired 2026-10-05).
- ~~Pluma–Speicher~~, ~~double-scaled SUSY SYK (BBNR)~~ (acquired 2026-10-05). Feng–Tian–Wei is still missing.
- New since 2026-10-05, not acquired:
  - ~~Penington–Witten, arXiv:2412.15549~~ (acquired 2026-10-27; `literature/notes/penington_witten_2024.md`).
  - ~~Lee–Lee, arXiv:1911.10391~~ (acquired 2026-10-27; `literature/notes/lee_lee_2019.md`). *Correction:* this entry
    said Lee–Lee supply "exchange data". They do not. They give multiplets, recombination rules and characters, and G. Lin
    cites them only for notation and channel sums.
  - **New gap:** the N=2 (su(1,1|1)) super-6j / exchange matrix Λ_Ξ for G. Lin's OTOC route (his eq. 6.8). G. Lin says
    it is not yet known; we have no source for it.
  - Aguilar-Gutierrez, arXiv:2510.11777: BPS Krylov complexity in N=2 DSSYK. Marginal for us (abstract checked
    only).
- Iniguez–Srednicki, arXiv:2305.15702, and Wang et al., arXiv:2310.20264: microcanonical truncations
  (Belin–Fu–La Rocca refs. [13], [11]).
- Feng–Tian–Wei (SYK spectrum q-Gaussian): Liu–Shen ref. [14]. (Pluma–Speicher, ref. [27], acquired 2026-10-05.)
- *Watch:* Khamnei–Papadodimas, unpublished, cited by Belin–Fu–La Rocca [53] for "few-fermion operators are not free
  after BPS projection in SUSY SYK". This overlaps qualitatively with our non-freeness results.
- ~~Double-scaled supersymmetric SYK~~ (acquired 2026-10-05 as BBNR). An earlier version of this entry gave the fourth
  author as "Torrents"; Berkooz–Mamroud ref. [66] has Raz.
- Saad–Shenker–Stanford (1903.11115) and Eynard–Orantin (math-ph/0702045), cited by the draft.
- Two citations in the draft (MP 2403.05241; DistinguishMicro 2108.00011) have unverified or possibly wrong metadata.

## 7. Where the project sits (2026-09-30)

The project combines three literatures that have not been combined before, as far as our collection and a brief
search (2026-09-30) show.

| Strand | Known in the literature | What the project adds (status) |
|---|---|---|
| **Fortuity / BPS chaos** (CV, CCSY, LMRS, Chen–Lin–Shenker, Belin–Fu–La Rocca, Miyahara–Shibuya, CCMP) | Projected operators diagnose BPS chaos. The Haar-orientation hypothesis gives GUE/Wachter statistics. In SUSY SYK, simple projected operators are not GUE (LMRS) and not free (Khamnei–Papadodimas, unpublished). | The uplift decoder recast as the principal-angle law of two correlated BPS spaces B, B′ = UB, with B′ the BPS space of the sign-flipped model (derived). The exact deviation from Haar orientation at second order is the U(N) size spectrum (derived + verified). SYK wall atoms exceed the Haar-generic count (num + Belin Thm. A.2). |
| **DSSYK chords / operator size** (BINT, Berkooz–Mamroud, Lin 2022, BLY) | Matter-chord crossing weights q^Δ. Chord number = size/q. Two-sided light probes measure length. N=2 super-chords and the zero-temperature two-point function. | The decoder variance as an exact finite-N two-sided light-probe correlator: r = 1 − 4⟨k(N+1−k)⟩/(N(N+1)), w₀ = a exactly (derived + verified). Size per chord unit = p − 1 in N=2 SYK (D3; derived for H, numerical for P). The Z₂ arc rule for a fixed parity probe (derived; standard per chord, a fixed-probe variant for m ≥ 4). |
| **Free probability / deviation** (Collins, Kunisky, Speicher, Liu–Shen, Belin–Fu–La Rocca) | Free compression = Wachter. Mixed q-Gaussians and ε-freeness for SYK Hamiltonians. Haar ⇒ asymptotic freeness. | Freeness = flat transmission (exact finite-N characterization for U(N)-invariant ensembles; derived + verified). The two-model fourth moment predicted from single-model data to ≤ 0.02 % (num). The exact decoder T₄ decomposition, with a four-point remainder R (derived + num; R open). |

**Honest assessment.**
- The chord rules themselves (per-chord x = matter weight, chord number ↔ size) are **standard**.
- The representation-theoretic layer is **plausibly new**: exact finite-N sum rules for the BPS projector, the
  transmission spectrum, freeness as flat transmission, and the T₄ decomposition.
- So is its application to the fortuity/uplift decoder.
- Closest competitors: Belin–Fu–La Rocca (the Haar side of the same question) and the unpublished
  Khamnei–Papadodimas work (non-freeness of projected simple operators in SUSY SYK).

## 8. Tools for the T₄ / R program (2026-10-05)

**Question.** What is the decoder fourth moment T₄ = τ(PP′PP′), P′ = UPU, U = (−1)^{n_i}? In particular, what is its
remainder R beyond the U(N)-twirl part? (`research/pdfs/progress_2026-09-29.pdf` §6–7; `docs/derivations.md` D5.)

| Tool | What it gives | Status / limits |
|---|---|---|
| Exact finite-N U(N) representation theory (project) | the size spectrum w_k, transmissions φ_k, sum rules, the T₄ = a² + 2a(T₂ − a) + S + R split | exact; leaves R as a genuine four-point quantity |
| OU-correlated couplings + second quantization (BKS Thm. 2.11; Pluma–Speicher Cor. 3.2; finite-N classical Mehler) | T₂^OU(x) = Σ κ_n x^{2n}; R^OU(x) = Σ ρ_n x^{2n} with R^OU(0) = 0 exactly; a one-parameter interpolation between independent and identical BPS spaces | T₄^OU is out of OU-time order, so the Markov kernel is not enough. Finite-N numerics (Phase A1) are needed |
| ℤ₂ four-arc chord rule (derived in project, D5) | the decoder and the OU four-arc function share one structure: adjacent-arc chords weigh x (x_s = 1 − 2p/N for the decoder), opposite-arc chords weigh 1 | the decoder's B-chords are structured (all contain site i), unlike OU's. The difference is what A2 measures |
| BINT crossed 4-pt / R-matrix | the closed-form crossing kernel for *random* matter | not directly applicable: random matter weighs opposite-arc chords q̃² rather than 1 |
| BBNR / BLY N=2 chord Hilbert space | the zero-energy HH state; the two-point function ⟨q^{Δ N̂}⟩ = T₂^OU at x = q^Δ; the one-particle wormhole for crossed functions | the crossed function is not evaluated in the literature. The finite-λ computation is open (our action item) |
| Super-JT (LMRS; G. Lin 2026) | zero-energy OTOC at equal Δ (exact, N=2); large-Δ suppression | the conformal regime overshoots our probe's 2-pt function by 15–60 % at N′ ≤ 16. Use as a limiting check only |
| Quantum geometric tensor (CCMP App. A) | small-angle OU: T₄ − 2T₂ + a ≈ (4ε²/D) E[Tr g_tot² + 2Σ Tr g_IJ²] | metric only. Berry curvature (TOC − OTOC à la G. Lin) needs two distinct directions in a loop, not isotropic noise |

**Recommendation.**
- Work at finite N with the OU family as the organizing interpolation, since it is exact and computable.
- Use the ℤ₂ chord rule to interpret the result. Compare against BBNR/BLY double scaling where it can be evaluated.
- Treat super-JT results as the triple-scaled limiting check, not as the quantitative model.
