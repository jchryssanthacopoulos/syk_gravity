# Literature synthesis

*Last revised 2026-09-28 (second pass, after adding 9 papers: PSSY, Turiaci–Witten, Boruch–Iliesiu–Yan, HITZ,
Heydeman–Turiaci–Zhao, FGMS, Y. Chen, Collins–Matsumoto–Novak, Berkooz–Mamroud). Per-paper notes are in
`literature/notes/`; BibTeX in `literature/bibliography.bib` (26 entries).*

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
                                                                                             fragility)
 ─── random matrices / free probability ────────────────────────────────────────────────────────▲
 Collins 2004 (exact Jacobi law) ─ Collins–Matsumoto–Novak 2021 (Weingarten 1/N = monotone walks) │
 Kunisky 2023 (MANOVA theorems, error identity) ─ Haikin–Zamir–Gavish 2017 ─ Dubbs–Edelman 2015    │
 Demni–Hamdi–Hmidi 2012 (free Jacobi/liberation) ─ Liu–Chen–Balents 2017 ─ Bouchard 2024 (TR) ────┘
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

## 6. Gaps in the literature collection

Still missing and worth adding:
- LMRS companion, arXiv:2207.00407 (*Holography for people with no time*).
- Boruch, Lin, Yan, *Exploring supersymmetric wormholes in N=2 SYK with chords*, JHEP 12 (2023) 151 — likely the
  most direct route to BPS-projected correlators of fixed-size operators at finite N.
- Double-scaled supersymmetric SYK (Berkooz, Brukner, Narovlansky, Torrents; ref. [66] of Berkooz–Mamroud).
- Saad–Shenker–Stanford (1903.11115) and Eynard–Orantin (math-ph/0702045), cited by the draft.
- Two citations in the draft (MP 2403.05241; DistinguishMicro 2108.00011) have unverified or possibly wrong metadata.
