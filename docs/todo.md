# To-do

*Last revised 2026-10-27 (paper roadmap checklist added at the top); 2026-10-26 (D10 uplift anywhere in the window; KW 2017 note); 2026-10-25 (D9 chords for the two-flavor uplift); 2026-10-23 (D8 two-flavor uplift); earlier 2026-10-12 (D6 uplift/localization), 2026-09-30 (literature pass), 2026-09-28 (post-archaeology). Linked questions refer to `docs/research_questions.md`.*
Priority: **P0** blocks everything / **P1** next / **P2** later.

## Paper roadmap (running checklist; started 2026-10-27)

**Goal.** A paper on the gravitational (chord) interpretation of fragility and fortuity in N=2 SYK. It has three parts:
1. the one-flavor uplift in every sector (fortuitous states only);
2. its chord picture;
3. the same for the two-flavor model (fortuitous and monotone states).
The payoff is a distinction between fortuitous and monotone states inside a bona fide gravity (chord) calculation.

### Done
- [x] One-flavor uplift at half filling: exact formula, decoder, dressing, walls, T₄ = T₂ − H (D6).
- [x] One-flavor uplift in every window sector (D10).
  - General theorem on admissible sources; three regimes (central / resonant / edge).
  - Exact ranks of CV's δ; exceptional states = radical, with the exchange.
  - Jump maps for q = 3, 5. Verified to N′ = 16.
- [x] One-flavor chords, central sectors, first order (D7): ASC measure, time integral, free compound Poisson law.
- [x] One-flavor chords, resonant sectors (D11):
  - layer 1: the heavy-tailed first order, K^{−1/2}, converging with N′;
  - layer 3: R = AA†, rank fraction e^{−10λ/9}, a free (Haar) model of the exact decoder, and the full chord model.
- [x] Two-flavor uplift (D8): obstructed mapping cones, protection identity [E, V] = K, counting laws, closed forms
  with self-energies.
- [x] Two-flavor chords at first order (D9): dictionary, window, emission, resonances at s′ = ±1/5, the tower does not
  dress (numerically).
- [x] Literature: Kanazawa–Wettig 2017 note; synthesis updated.

### A. One-flavor chords: finish everywhere
- [ ] A1. The exact (self-dressed) decoder of the central sectors from chords: R and the h–g mixing block (0.008–0.019
  of T₂). Same machinery as D11 layer 3.
- [ ] A2. Layer 2, the finite-N hard edge at the window edge: residual tail deficit (12–26 %), extreme modes, KW
  GOE/GSE class. Appendix level.
- [ ] A3. A chord derivation of R's spectrum from the continuum→continuum N=2 matter vertex. It replaces the c_R/E_w
  closure.
- [ ] A4. The deep-tail alignment of range R with the softest edge modes (the free model is up to 45 % low for
  λ < 0.03 at N′ = 16). Caveat plus appendix.

### B. Two-flavor: bring to parity with one-flavor
- [ ] B1. D10-style regime classification in the 2D window (gapped / resonant / obstructed), with ranks and a jump map.
- [ ] B2. Resonant layers 1 and 3 for two flavors at s′ = ±1/5. This includes the chord computation of the
  self-energies 𝒮_G, the transfer 𝒵 and the leaks α, γ (open in D8/D9).
- [ ] B3. Larger two-flavor sizes, e.g. N′ = 10 at resonant parity, via sparse or partial methods.

### C. The fortuitous vs monotone distinction (the crux)
- [ ] C1. Derive the tower's protection inside the chord theory: resum the twisted chords (V-descendants) into a chord
  Hilbert space in which emission into the continuum vanishes. Highest priority and highest risk.
- [ ] C2. A quantitative side-by-side in chord variables. Fortuitous: dressing distribution, resonances, decoder,
  fragility. Monotone: zero dressing, slot weights (N+1−n)/(N+1), pair absorption c = n(N′−n)/N′². One table, one figure.
- [ ] C3. The Schwarzian (triple-scaled) limit of both behaviours; connect to super-JT (HITZ, LMRS).
- [ ] C4. Pressure-test the interpretation: is "tower = graviton gas" justified (check CCSY, HTZ)? Do any one-flavor
  states, e.g. the exceptional states, behave tower-like?

### D. Paper assembly
- [ ] D1. Unified notation and dictionary across models: charges s, matter δ, regimes, and the admissible-source
  theorem stated once.
- [ ] D2. Consolidated figures. Main text: jump maps (both models), regime diagram, dressing vs protection.
- [ ] D3. Acquire the literature: Penington–Witten (2412.15549), Lee–Lee (1911.10391), LMRS companion (2207.00407).
- [ ] D4. Optional: the pinned pre-registered N′ = 17 chord test (D7). Run only when asked.
- [ ] D5. Draft the paper. Half filling as the worked example; "the whole window" as one main-text section; details in
  appendices.

**Suggested order.** C1 → B1, B2 → A1 → C2, C3 → A2, A3 → D.

---

## A. Environment & reproducibility (P0)

- [ ] **Fix the Python environment.** The Anaconda base env has NumPy 2.2.4 with SciPy 1.13.1 / matplotlib compiled
      against NumPy 1.x, so `import scipy.sparse` and `import matplotlib` fail. Every script except the pure-NumPy
      checks is currently unrunnable there. A clean venv (numpy 2.5.3, scipy 1.18.1) runs `decoder_fast.py` and
      `protected_atoms.py` fine. Create a project venv and pin versions in `requirements.txt`.
- [ ] `requirements.txt` lists `cvxpy` and `pypdf`, which no code imports. Remove or justify; add `sympy` (audit scripts).
- [ ] Initialize git (the directory is not a repository) and make a first commit of the current state *before* any refactor.
- [ ] Record the code version (commit hash) in every future results file.

## B. Data hygiene (P0)

- [ ] `results/data/softtail_scaling.csv` mixes two schemas (12 vs 14 columns) and duplicates the N=15, seed 0
      row. Split into versioned files or add the missing columns; do not overwrite. Document which script version
      wrote which rows.
- [ ] The N = 12, 13 soft-tail numbers quoted in the draft (9%, 9%) are not in `results/`. Regenerate and store them.
- [ ] `results/data/protected_atoms.csv` has only N = 15, 16. Add N = 8…14 (the N = 8–11 values from the audit
      smoke test: 45/81, 90/162, 80/243, 160/486).
- [ ] Add a `results/README.md` recording, for each file: script, parameters, seeds, date.

## C. Code audit follow-ups (P1)

- [ ] Consolidate the **five** independent exterior-algebra implementations (`reproduce_table1.py`, `decoder_fast.py`,
      `protected_atoms.py`, `q_scan.py`, `KT_histogram.py`) into one module (e.g. `src/syk_forms.py`: subsets,
      signs, wedge matrices, harmonic basis). Add a cross-check test that all give identical spectra for a fixed seed.
- [ ] `src/KT_histogram.py` computes B_old† P_{B_new} B_old (compression of K_T onto old BPS states), **not** the
      O_T spectrum of the draft. Its Wachter overlay uses a = d_old/C(N,p), b = (N+1−p)/(N+1), which is not a derived
      null model for that operator. Label it clearly, and derive the correct null model (Q7) or remove the overlay.
      Figures `results/figures/N9_p4_down.pdf`, `N9_p5_down.pdf`, `my_spectrum.pdf` come from it.
- [ ] `decoder_chaos.py` docstring refers to "Figure 12 / Section 8.4 / Eq. 98–100" of an older numbering. CV now
      has it as Fig. 13 / eq. 105 with 500 realizations. Update.
- [ ] `decoder_fast.tail_stats`: κ is computed on λ ≥ 1e−9 (atoms dropped silently). Report the atom count
      alongside κ, and use a threshold-independent tail measure (Q5b).
- [ ] Label GUE ⟨r⟩ benchmarks consistently (surmise 0.60266 vs large-N numerics 0.5996).
- [ ] Write the missing generators for the draft's figures, or mark those figures as unreproducible:
      density/CDF vs Wachter (N = 6…10), soft_tail, sourcing, sff_ramp, moment_comparison_N10, rhodec_vs_syk,
      softtail_scan, chaotic_vs_nonchaotic. `research/tex/report_figs/` is empty.
- [ ] Add `tests/`: particle–hole λ ↔ 1−λ at even N′ (per realization); m₁ = 1/2 at even N′ half filling;
      disorder-averaged m₁ = b; d = dim B^P_N + dim B^{P−1}_N; atom rank bounds.

## D. Verification of existing derivations (P1)

- [x] Wachter moments m₁–m₃, S-transform, spectral curve (hand + sympy).
- [x] Weingarten m_k^{(1)}, k ≤ 5, and the W₁ closed form (`scripts/audit/weingarten_W1_check.py`).
- [x] Kunisky/Chebyshev kernel identity; eq. (errdef) typo (`scripts/audit/kernel_identity_check.py`).
- [ ] Derive the W₁ normalization analytically (expected (D/d)² from 1/d² vs 1/D² expansion, up to sign).
- [ ] Derive and numerically verify the constrained generalized eigenproblem replacing the false self-energy
      formula O_T = (1+C′†GC′)^{−1} (Q5a).
- [ ] Verify the bordered-Toeplitz claims (Lanczos coefficients from exact spectra; "length two"; P₄). Expected
      conclusion: they cannot capture a continuous sub-edge tail.
- [ ] Verify the cylinder W_{0,2} normalization and compute the SFF of O_T in SYK (not a Jacobi surrogate).

## E. Research (P1 → P2)

- [ ] **Q1 (P1):** LMRS zero-energy 2-point function of n_{N′} − ⟨n⟩ in the fixed-P BPS sector; compare to
      Var(λ)(N) for N′ = 8…16 (exact) and beyond (matrix-free).
- [ ] **Q2 (P1):** principal angles between B and B⁰ vs N (more seeds, N′ up to 15–16); first-order perturbation
      theory in the new couplings; free-Jacobi-process fit.
- [ ] **Q3 (P1):** explain the atom excess over the rank bound (N = 10–13) via constrained rooms.
- [ ] **Q6 (P2):** Thouless time of O_T; decoder along the Miyahara–Shibuya Q_g interpolation; two-flavor model.
- [ ] **Q7 (P2):** reconcile the three decoder spectra (CV state-dependent, uniform O_T, compression on B_old).
- [ ] **Q8/Q9 (P2):** q-scan once Q1 is settled; single-mode entanglement framing.
- [x] **Q10 (P1):** size ↔ chord dictionary derived (D4, `research/notes/chord_size_dictionary_2026-10-02.md`).
- [ ] **Q10 follow-ups (P2):** finite-N HH recursion for πₙ (n ≥ 2); smearing law; p = 5 exactness of P₍₁₎ ∝ :H:;
      j_R ≠ 0 joint distribution.
- [x] **Q11 (P1):** decoder vs OU vs free liberation; R decomposed (D5, `research/notes/decoder_T4_fourarc_2026-10-05.md`).
- [x] **Q2/Q11 (P1), 2026-10-12:** exact uplift, orthogonal splitting and closed-form block decoder (D6 §§1–5);
      localization = heterogeneity deficit, predicted from the old model's first-order uplift (D6 §6,
      `research/notes/uplift_localization_2026-10-12.md`). Wall hypothesis refuted at N′ = 15.
- [x] **D6 follow-up, 2026-10-14:** chord computation of the first-order dressing spectrum (D7, report
      `research/pdfs/progress_2026-10-14.pdf`, data `results/data/chord_dressing_2026-10-14/`).
- [ ] **D7 follow-ups (P1):**
  - run the pre-registered N′ = 17 test (`prereg_N17.json`; `scripts/run_chord_dressing.py --Np 17 --seeds 9000 9001
    --neta 4 --draws 2 --no-third-moment` under memwatch, then `scripts/analyze_bps_heterogeneity.py --Np 17`);
  - complete the stopped N′ = 16 runs (second realization of the decomposition; heterogeneity analysis);
  - large-N fate of the coherent BPS heterogeneity; relate B_c (single-fermion BPS weight) to the decoder's
    localization;
  - a model of the soft lower edge (operator-valued free probability with the exact covariance of Σ η_I A_I);
  - N=2 matter–matter crossing (R-matrix); low priority (crossed diagrams ≤ 0.5 %).
- [ ] **D6 follow-ups (P1):**
  - why first order preserves H = ⟨y(1−y)⟩ (approximate cancellation; no symmetry found);
  - the odd-N′ residual of H (−0.001, −0.003, −0.004 at N′ = 11, 13, 15): first-order N′ = 17 (old model N = 16; cheap)
    plus one exact N′ = 17 point if memory allows;
  - the deterministic atom excess over the refined bound (+1, +2, +14, +28 at N′ = 11–14).
- [x] **Q12 Stage 1, 2026-10-23:** exact two-flavor pair uplift (D8; report `research/pdfs/progress_2026-10-23.pdf`;
      data `results/data/two_flavor_2026-10-23/`).
- [ ] **Q12 / D8 follow-ups:**
  - ~~**Stage 2 (P1):** chord computation of the E (ψ†χ), K (ψ†ψ†), Z (ψ†) two-point functions and their interference;
    first check the literature for two-flavor N=2 chord rules~~ → done 2026-10-25 (D9; report
    `research/pdfs/progress_2026-10-25.pdf`; data `results/data/two_flavor_chords_2026-10-25/`): two-flavor chord dictionary,
    black-hole window with gapless edges, emission resonances, tower = no emission (twisted chord = V insertion);
  - **Stage 2 follow-ups (P1):** chord evaluation of 𝒵, α, γ (exact uplift from chords); a chord Hilbert space for graviton
    states (resum twisted chords); per-sector chord prediction of the mean exact dressing (slot weights + interference);
    N' = 10 at resonant parity with sparse methods (predicted: the J-resonant bulk sectors (3,2), (4,4), (5,6) and their
    mirrors (4,2), (5,4), (6,6) far above the s' = 0 sectors (4,3), (5,5), (6,7), with the tower bare in (4,4), (6,6));
  - ~~proofs: four-slot injectivity; closed form of d₂~~ → done 2026-10-25 (D8 §9): d₂ closed form via homological
    perturbation; injectivity is generic but fails at special couplings (pure-dressing states);
  - proofs: the room laws L1/L2 (generic position in the refined rooms); the geometry of the pure-dressing set;
  - the N = 5 coincidences (three one-dimensional excesses in near-central sectors, both models);
  - ~~closed-form CV slot decoders~~ → pencils in D8 §9.6; still open: simplifying them (push-through) and the
    dressing-cost spectra κ_T of fortuitous states vs CV Fig. 8;
  - N = 8 → 9 needs sparse/partial methods (largest sector 15876).
- [ ] **Q11 follow-ups (P1):**
  - ~~analytic model of the localization term~~ → D6 §6 (semi-analytic: old-model first-order uplift);
  - ~~N′ = 15–16 for its large-N trend~~ → N′ = 15: L flat (D6 §6);
  - two-colour chord computation of the isotropic (OU) four-arc function in BBNR/BLY double scaling (the
    mismatched-wormhole overlap);
  - Krawtchouk inversion of the parity family (u_t grading, D5 §1).

- [ ] **Q13 follow-ups (P1, D10):**
  - the N′ = 9 radical excess: test whether rank B_η on H(8, −2) is 24 for every η, including structured choices
    (rank-deficient η);
  - large-N resonant tail: extend N′ = 16 (done 2026-10-26: α = 0.48 ± 0.03) to N′ = 18 with sparse/partial methods
    (only the lowest edge modes matter; rank decisions must use the known dimensions, D10 §5 caution);
  - chord amplitude of the tail beyond the rank-one approximation (free compound Poisson with heavy-tailed jumps);
  - an analytic model for the exact resonant decoder edge λ_min ≈ 10⁻² (the spectral measure of R along the
    resonant modes);
  - general q: exceptional counts beyond parity (q = 5, N′ = 19, 23 need sparse methods).

## F. Writing (P2)

- [x] Progress report 2026-10-26 (`research/pdfs/progress_2026-10-26.pdf`, source `research/tex/progress_2026-10-26.tex`):
      the uplift anywhere in the window (D10); jump map, exact ranks and exceptional radicals, edge resonance and its
      chord tail; KW 2017 note.
- [x] Progress report 2026-10-23 (`research/pdfs/progress_2026-10-23.pdf`, source `research/tex/progress_2026-10-23.tex`):
      two-flavor uplift; obstruction theorem, protection identity, tower bare, CV Fig. 6 derived from counting laws,
      bare-state count, first-order interfering two-point functions (D8).
- [x] Progress report 2026-10-14 (`research/pdfs/progress_2026-10-14.pdf`, source `research/tex/progress_2026-10-14.tex`):
      the dressing spectrum from chords (closed-form BPS→continuum two-point function, ₂φ₁ time integral, free compound
      Poisson law, error budget, sum rule and BPS heterogeneity, pre-registered N′ = 17 test) (D7).
- [x] Progress report 2026-10-13 (`research/pdfs/progress_2026-10-13.pdf`, source `research/tex/progress_2026-10-13.tex`):
      exact uplift, orthogonal splitting, closed-form decoder, dressing cost as integrated two-point function, walls,
      localization as heterogeneity deficit with the first-order prediction and the N′ = 15/16 pre-registered tests (D6).
- [x] Progress report 2026-10-05 (`research/pdfs/progress_2026-10-05.pdf`, source `research/tex/progress_2026-10-05.tex`):
      chord ↔ size dictionary (D4) and the decoder fourth moment as liberation + localization (D5).

- [ ] Rewrite the draft following the "safe core" (see `research/notes/gravitational_dual_fragility_audit.md` §11):
      exact one-flavor geometry → numerics → Haar null model (clearly labelled) → LMRS calculation (Q1) → only then
      any gravity interpretation. Move the two-flavor material to a separate control section.
- [ ] Fix the bibliography: Chang–Chen–Sia–Yang for 2412.06902; add authors to 2608.12160; verify MP,
      DistinguishMicro, N2matrix metadata.
- [ ] Keep `docs/derivations.md` for the verified derivations (block BPS equations, norm identities, atom
      characterization and rank bounds, Weingarten/W₁ for the null model).

## G. Literature (P2)

- [x] Added notes (2026-09-28) for FGMS, Turiaci–Witten, HTZ, Boruch–Iliesiu–Yan, PSSY, HITZ, Y. Chen,
      Collins–Matsumoto–Novak (Weingarten), Berkooz–Mamroud (DSSYK chords).
- [x] Boruch–Lin–Yan (arXiv:2308.16283) acquired and noted 2026-09-28 (`literature/notes/boruch_lin_yan_2023.md`).
- [x] 2026-10-26: Kanazawa–Wettig 2017 acquired and noted (`literature/notes/kanazawa_wettig_2017.md`); synthesis §3 item 13.
- [ ] Still to acquire: LMRS companion (2207.00407); double-scaled supersymmetric SYK (Berkooz et al.); SSS; Eynard–Orantin.
- [x] 2026-09-30: notes for Belin–Fu–La Rocca 2026, Lin 2022, Liu–Shen 2026; synthesis §7 (positioning).
- [ ] Acquire BINT 2018 (1811.02584, original matter-chord rule), Iniguez–Srednicki (2305.15702), Wang et al.
      (2310.20264), Pluma–Speicher, Feng–Tian–Wei.
- [ ] Watch for Khamnei–Papadodimas (non-freeness of projected simple operators in SUSY SYK; unpublished, cited by
      Belin–Fu–La Rocca [53]).
- [ ] Next report: cite Berkooz–Mamroud eq. 2.28 / BINT for the per-chord weight (progress report §3.2), Lin eq. 53/57
      for size ↔ length, Liu–Shen for inter-family crossings, and Belin–Fu–La Rocca for the Haar hypothesis.
