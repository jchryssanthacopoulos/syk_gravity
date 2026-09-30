# Literature synthesis

*Last revised 2026-09-30 (third pass: added Belin–Fu–La Rocca 2026, Lin 2022, Liu–Shen 2026; new §7 on where the
project sits). Second pass 2026-09-28 added 9 papers: PSSY, Turiaci–Witten, Boruch–Iliesiu–Yan, HITZ,
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
- LMRS companion, arXiv:2207.00407 (*Holography for people with no time*).
- ~~Boruch, Lin, Yan~~ (acquired 2026-09-28).
- Berkooz, Isachenkov, Narovlansky, Torrents, *Towards a full solution of the large N double-scaled SYK model*,
  arXiv:1811.02584: the original matter-chord crossing rule q̃ = e^{−2pp̃/N}.
- Iniguez–Srednicki, arXiv:2305.15702, and Wang et al., arXiv:2310.20264: microcanonical truncations
  (Belin–Fu–La Rocca refs. [13], [11]).
- Pluma–Speicher (independent SYK → q-Gaussian systems) and Feng–Tian–Wei (SYK spectrum q-Gaussian): Liu–Shen refs.
  [27], [14].
- *Watch:* Khamnei–Papadodimas, unpublished, cited by Belin–Fu–La Rocca [53] for "few-fermion operators are not free
  after BPS projection in SUSY SYK". This overlaps qualitatively with our non-freeness results.
- Double-scaled supersymmetric SYK (Berkooz, Brukner, Narovlansky, Torrents; ref. [66] of Berkooz–Mamroud).
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
