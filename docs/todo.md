# To-do

*Last revised 2026-10-23 (D8 two-flavor uplift); earlier 2026-10-12 (D6 uplift/localization), 2026-09-30 (literature pass), 2026-09-28 (post-archaeology). Linked questions refer to `docs/research_questions.md`.*
Priority: **P0** blocks everything / **P1** next / **P2** later.

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
  - **Stage 2 (P1):** chord computation of the E (ψ†χ), K (ψ†ψ†), Z (ψ†) two-point functions and their interference;
    first check the literature for two-flavor N=2 chord rules;
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

## F. Writing (P2)

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
- [ ] Still to acquire: LMRS companion (2207.00407); double-scaled supersymmetric SYK (Berkooz et al.); SSS; Eynard–Orantin.
- [x] 2026-09-30: notes for Belin–Fu–La Rocca 2026, Lin 2022, Liu–Shen 2026; synthesis §7 (positioning).
- [ ] Acquire BINT 2018 (1811.02584, original matter-chord rule), Iniguez–Srednicki (2305.15702), Wang et al.
      (2310.20264), Pluma–Speicher, Feng–Tian–Wei.
- [ ] Watch for Khamnei–Papadodimas (non-freeness of projected simple operators in SUSY SYK; unpublished, cited by
      Belin–Fu–La Rocca [53]).
- [ ] Next report: cite Berkooz–Mamroud eq. 2.28 / BINT for the per-chord weight (progress report §3.2), Lin eq. 53/57
      for size ↔ length, Liu–Shen for inter-family crossings, and Belin–Fu–La Rocca for the Haar hypothesis.
