# Project archaeology report

*2026-09-28. Scope: every file in the repository. All 17 PDFs in `papers/` were read (large ones partly
skimmed, noted per paper). The current draft was read from its TeX source, the prior audit PDF in full, and all
10 scripts and the results data. Three independent checks were run (`scripts/audit/`).*

Companion documents:
- per-paper notes: `literature/notes/*.md`; synthesis and map: `literature/synthesis.md`; BibTeX: `literature/bibliography.bib`
- claim-by-claim audit of the current draft: `research/notes/gravitational_dual_fragility_audit.md`
- `docs/research_questions.md`, `docs/todo.md`

---

## 0. Inventory

| Location | Contents | Provenance |
|---|---|---|
| `CLAUDE.md`, `README.md` | Project instructions; one-line README | user |
| `papers/` (17 PDFs) | Primary literature, incl. the project's own paper (Chryssanthacopoulos–Vegh 2026) | external |
| `research/tex/gravitational_dual_fragility.tex` + `research/pdfs/…pdf` | **Current draft** "Fragile Fortuity: The Gravity of BPS Overlaps" (30 pp, no authors) | Claude-generated |
| `research/pdfs/decoder_notes_audit.pdf` | 9-pp "Technical audit of the decoder–gravity notes" | Claude-generated (Sept 2026) |
| `research/tex/report_figs/` | empty (figures referenced by the draft are missing) | — |
| `src/` (10 scripts) | decoder ED, atoms, q-scans, CV Table 1 / Fig. 13 reproductions, Johnson Figs. 4–5 | mostly Claude-generated; `KT_histogram.py` apparently user-written ("WIP convention") |
| `results/data/` | `softtail_scaling.csv`, `protected_atoms.csv`, spectra N=14,15 (3 seeds each) | from `decoder_fast.py`, `protected_atoms.py` |
| `results/figures/` | Wachter comparisons (N 12→13, 13→14, 14→15), KT histograms (N=9), Johnson repros | from the scripts above |
| `docs/`, `literature/` | were empty (0-byte placeholders); now populated by this audit | — |
| `scripts/`, `tests/` | were empty; `scripts/audit/` added by this audit | — |

Not a git repository. No `research/notes/` or `research/synthesis.md` existed before this audit.

## 1. Research question

The material, read together, is about **one operator**: the down-channel uplift decoder of one-flavor N=2 SYK,

  O_T = P_B Π_T P_B,  Π_T = 1 − n_{N′},  P_B = projector onto BPS states of the (N+1)-fermion theory at charge P.

Its eigenvalues λ ∈ [0,1] are squared cosines of principal angles between the enlarged BPS space and the "slot"
(new mode empty). CV introduced it to quantify **metric fragility**: how well BPS states survive N → N+1. They showed
its interior level statistics are GUE (**BPS chaos** in the LMRS sense). They left a gravitational construction as
future work.

The project's aim (CLAUDE.md, current draft) is to understand the **spectral density** of O_T:
1. identify the right random-matrix law and quantify the deviation from it;
2. derive it from gravity (N=2 super-JT, chord diagrams, non-crossing partitions);
3. explain what the "free" description misses: atoms at λ = 0, 1 and a macroscopic "soft tail".

After the audit, the sharpest version of the question is: **is the decoder spectrum governed by a Haar/free-projector
law (Wachter), or by the LMRS projected-operator law of a light neutral bilinear, and what is the structure of the
deviations?** (See `docs/research_questions.md` Q1–Q2.)

## 2. Current theoretical understanding

### 2a. Established in the literature
- One-flavor N=2 SYK as a wedge with a generic 3-form. BPS = harmonic forms ≅ Q-cohomology. R-charge concentration
  in the fortuity window; all BPS states fortuitous (CCSY 2024; CV §2). Central BPS dimensions 3^{(N−1)/2} (N odd)
  and 2·3^{N/2−1}, 3^{N/2−1} (N even) (Berry Table 2; CV Table 1).
- Long exact sequence relating N and N+1; at least one connecting map vanishes in the window; uplift along walks;
  d = dim B^P_N + dim B^{P−1}_N when connecting maps vanish (CV eqs. 24–27, §4).
- Decoder D_T, K_T, O_T, and state-dependent fidelity / dressing cost (CV eqs. 42–63, 101–102). Exact-lift loci
  obey generic-position counting inside constrained rooms (CV Table 1). ⟨r⟩ = 0.599 ± 0.001 (GUE) at N′=12 (CV eq. 105).
- Wachter law = Bern(a) ⊠ Bern(b); exact finite-size Jacobi law for Haar P (Collins Thm. 2.2). Freeness criterion
  and exact error identity for products of projections (Kunisky Thm. 1.5, Lemma 4.6). Bordered-Toeplitz Jacobi data
  (Dubbs–Edelman).
- LMRS: zero-energy correlators of Ô = POP in N=2 super-JT depend on Δ and R-charge. Crossings on the disk are
  suppressed only at large Δ. e^{S0} counts BPS states (LMRS eqs. 65–66, 82, 159, 169; Berry §6).

### 2b. Derived in the project material (and verified in this audit)
- Block BPS equations for Ψ = a + b∧e (draft eq. bpsfour; = CV eqs. 20–22, 96).
- Characterization of exact atoms: λ=1 ⇔ a ∈ B^P_N ∩ ker(D∧); λ=0 ⇔ b ∈ B^{P−1}_N ∩ ker ι_D̄.
- Exact norm identities λ = E_b/(E_b+σ_a²), 1−λ = E_a/(E_a+σ_b²) (true but near-tautological).
- For the **Haar null model**: exact Weingarten moments, genus-one coefficients m_k^{(1)} and the closed form
  W₁(x) = −b(1−a)(1−b)x(x−1)/[(x−λ₋)(x−λ₊)]^{5/2} (verified k ≤ 5 here, k ≤ 6 in the prior audit).
- The Kunisky/Chebyshev kernel K_{k,ℓ} (an exact algebraic identity for any projector pair).
- Particle–hole pairing λ ↔ 1−λ at even N′, half filling (prior audit; consistent with m₁ = ½ exactly here).

### 2c. Numerical observations
- Interior spectrum of O_T is **not** Wachter. Weight is removed from mid-band (~1.7–2× depleted at N = 12–14) and
  piled up at both walls. KS to the (renormalized) Wachter interior is 0.12 (12→13) and 0.17 (14→15) over the full
  interior (`results/figures/decoder_wachter*.png`). Restricted to [λ₋, λ₊] it is 0.07–0.08
  (`softtail_scaling.csv`); that restriction hides the wall pile-up.
- Sub-edge population ~15% at N = 14–15 (repo data), ~10% at N = 12–13 (draft; reproduced by prior audit).
- δm₂ = m₂ − m₂^W ≈ 0.02–0.045 and not decreasing for N ≤ 12. Var(λ) exceeds ab(1−b).
- Atom fractions 55.6% (N = 8, 9) → 32.9% (10, 11) → 5.6% (12, 13) → 0 (14–16); counts ≈ rank lower bounds.
- Edge eigenvectors lie ≈ 97–99.6% in the LES reference subspace B⁰ (exploratory, N′ ≤ 12).
- GUE level statistics for the interior (CV). Tail ⟨r⟩ ≈ GUE (draft; no code).

### 2d. Hypotheses (plausible, untested)
- **H1:** the physical decoder variance decays as a power of N (LMRS, Δ = 1/3), unlike the exponentially small
  Wachter variance.
- **H2:** O_T is a partial liberation of B⁰ (whose spectrum is Bernoulli on {0,1}) toward Wachter; walls are the
  un-liberated remainder.
- **H3:** atom excess over rank bounds is explained by generic position inside constrained rooms.
- **CV conjecture:** metric fragility is a signature of BPS chaos.

### 2e. Speculative / insufficiently justified (from the draft)
Two-species EOW-brane model with e^{S0} = D; "BPS projection trivializes the Schwarzian"; "the super-JT action
reproduces the handles to all orders"; "crossing = handle"; self-energy form of O_T; "the tail is the overlap shadow
of the near-extremal throat"; length-two bordered-Toeplitz closed form; "the bulk disappears at large N";
"fragility = entanglement deficit". Details in the draft audit, §§1–9.

## 3. Literature map

See `literature/synthesis.md` §2 for the diagram and table. In brief:

- **Chryssanthacopoulos–Vegh 2026** — defines the decoder, diagnostics, LES/walk framework, GUE result → the whole project.
- **LMRS 2022** — the gravity calculation for BPS-projected operators; contradicts the draft's gravity model;
  predicts power-law N-dependence → Q1, Q4.
- **Chen–Lin–Shenker 2024** — LMRS criterion, strong/weak chaos, random-hyperplane argument (valid only for ν ≪ L)
  → interpretation, Q6.
- **Chen–Colin-Ellerin–Mamroud–Papadodimas 2026** — chord diagrams at Δ→∞, free commutator, SYK BPS numerics, UV
  imprint on Λ̂ → draft's chord section, Q1.
- **Chang–Chen–Sia–Yang 2024** — model, R-charge concentration, two-flavor monotone tower, LMRS test → model/control.
  (Mis-cited in the draft.)
- **Johnson 2026** — energy-side Wishart/Bessel model; repro scripts → contrast only.
- **Miyahara–Shibuya 2026** — BPS-subspace chaos↔integrability crossover with fortuitous states → test of CV conjecture (Q6).
- **Kunisky 2023** — freeness criterion; exact error identity → draft Sec. 7 (correctly cited, over-interpreted);
  open problem on smallest eigenvalue = soft tail.
- **Collins 2004** — exact Jacobi law; edge universality; large deviations → Haar null model.
- **Dubbs–Edelman 2015** — bordered-Toeplitz Jacobi data → draft's ρ_dec (which cannot hold a continuous tail).
- **Haikin–Zamir–Gavish 2017** — structured subspaces can be MANOVA → caution.
- **Demni–Hamdi–Hmidi 2012** — free Jacobi / liberation process → H2.
- **Liu–Chen–Balents 2017** — Wachter in SYK₂ entanglement → draft's EE bridge (parameter pairing reversed).
- **Bouchard 2024** — TR conventions → draft App. A.
- **Chang–Lin 2024 / 2022, Chang–Yin 2013** — fortuity background.

## 4. Existing derivations

| Derivation | Where | Assumptions | Correct? | Needs |
|---|---|---|---|---|
| Block BPS equations | draft §1.2; CV §3.1, §7.1 | one added mode, q=3 | ✅ (signs convention-dependent) | — |
| Hodge projector P_B = 1 − QGQ† − Q†GQ | draft eq. PB | — | ✅ | — |
| Self-energy form O_T = (1+C′†GC′)^{−1} | draft eq. KTself | b has no harmonic part; domain ker Q_N | ❌ | Replace by the constrained generalized eigenproblem (prior audit's proposal, unverified) |
| λ = E_b/(E_b+σ_a²), mirror identity | draft eqs. identity, identityup | eigenvector of O_T, b ∈ ker Q†, a ∈ ker Q | ✅ (tautological) | Stop calling the linear law a theorem |
| Tail = image of ρ_H | draft eq. tailimage | E_b ≈ eigenvalue of H_N; σ_a² constant | ❌ | Superseded by near-harmonic / B⁰ picture |
| Wachter law, moments, S-transform, curve | draft §1.4 | free projectors | ✅ | — |
| m₁ = b "exactly" | draft §1.4 | — | ⚠️ disorder-average only | — |
| Weingarten genus expansion, m_k^{(1)} | draft §2 | Haar P_B | ✅ (k ≤ 5 here) | Fix a=b=½ numeric |
| W₁ closed form via TR | draft §3.6, App. A | Haar; TR on Wachter curve | ✅ (Haar) | Derive normalization; drop "whole tower" |
| Amplitude dictionary from super-JT action | draft §4, App. B | e^{S0}=D; EOW slot branes; trivial Schwarzian | ❌ | Would need an LMRS-type computation (Q4) |
| Kunisky error recursion, K_{k,ℓ} | draft §5 (errdef–Kclosed) | any projectors | ✅ identity; eq. errdef typo | Reframe as bookkeeping |
| Bordered-Toeplitz ρ_dec, P₄ | draft eqs. borderedtoeplitz–P4 | "length two" | ❓/❌ | Cannot represent a continuous tail |
| d(N) closed form | draft eq. dclosed | N = old size | ✅ | Clarify N vs N′ |
| a(N) asymptotics | draft eq. aasy | — | ⚠️ even N only | parity factor 2/√3 |
| Band collapse λ± ≈ ½ ± √a | draft eq. collapse | Haar | ✅ | Interpretation wrong |
| Walk counting of Hodge sectors | CV App. A | Conj. A.1 | ✅ (conditional) | — |
| Asymptotic obstruction (codim>0 for odd N ≥ 15) | CV eq. 100 | Conj. A.1, generic position | conditional | — |

`docs/derivations.md` is still empty. The verified items above are the ones to transcribe there.

## 5. Existing numerical work

| Code | What it does | Status |
|---|---|---|
| `src/reproduce_table1.py` | CV Table 1 (exact-lift loci) | Runs (per structure); not re-run here |
| `src/decoder_chaos.py` | CV Fig. 13 / ⟨r⟩ | Imports `reproduce_table1`; docstring numbering outdated |
| `src/decoder_fast.py` | Sparse BPS basis, O_T spectrum via SVD of empty block, tail stats → `softtail_scaling.csv` | **Runs** (clean venv), checked N′ = 9, 10 |
| `src/decoder_wachter_check.py` | O_T vs Wachter histograms → `decoder_wachter*.png` | Needs matplotlib |
| `src/protected_atoms.py` | Exact atom counts via kernels → `protected_atoms.csv` | **Runs**; N = 8–11 checked here |
| `src/q_scan.py`, `q_scan_mf.py` | δm_k vs q; matrix-free Hutchinson estimator | Not run; no results stored |
| `src/KT_histogram.py` | **Different operator** (B_old† P_new B_old) with an underived Wachter overlay → `N9_p*_down.pdf`, `my_spectrum.pdf` | Needs clarification (Q7) |
| `src/reproduce_johnson_fig{4,5}.py` | Wishart repro of Johnson Figs. 4–5 | Consistent with Johnson eqs. 26, 36, 37 |

**Established numerically:** GUE interior statistics (CV); non-Wachter interior with wall pile-up; sub-edge
population ~10–15% at N = 12–15; atoms vanish for N ≥ 14 (N = 14–16, 1–2 seeds); d = dim B^P_N + dim B^{P−1}_N
(all N in the data); particle–hole pairing; exact norm identities.

**Untested or unreproducible from repo:** all figures in the draft (no generating code or data); SFF of the SYK
decoder (only a Jacobi surrogate was used); tail ⟨r⟩; Lanczos "length two"; σ_a² statistics; N = 12, 13 tail data;
q-dependence; anything beyond N′ = 16. Statistics are thin: 1–3 seeds per N for tail data.

**Environment:** the Anaconda base env is broken for scipy/matplotlib (NumPy ABI mismatch). Use a venv.

## 6. Open questions

Prioritized list with concrete tasks in `docs/research_questions.md`:
1. **Q1:** Haar (exponential) vs LMRS (power-law, Δ=1/3) scaling of the decoder variance. Compute the LMRS prediction.
2. **Q2:** partial-liberation / reference-subspace structure of O_T.
3. **Q3:** exact atom counts via generic position.
4. **Q4:** a legitimate gravitational description (e^{S0} ∝ d; Π_T as operator insertion; finite-Δ crossings).
5. **Q5:** correct eigenproblem and large-N fate of the soft tail, with a threshold-independent measure.
6. **Q6:** strength of BPS chaos (Thouless time) and the fragility↔chaos conjecture (integrable interpolation).
7. **Q7:** which decoder spectrum is the physical object (three different ones are in play).

## 7. Contradictions and uncertainties

**Inconsistent notation / conventions**
- N means the *old* size in the draft (d(12) = 729 = dim B^6_{13}) but the *target* size in the prior audit. This
  caused the audit's false "factor 3/2 error" claim. Adopt N (old), N′ = N+1 (target) everywhere.
- S₀ means ln D (§§2–4) and the BPS entropy ½ ln 3·N (eq. index) in the same draft.
- F_T, κ_T: state-dependent over old states (CV) vs uniform over new BPS space (draft).
- "Fragility": metric fragility (CV) vs fragility of monotone states under coupling perturbations (CCSY).
- Wachter parametrizations: (α,β) ∈ (0,1) (Kunisky, Collins) vs (a,b) ≥ 1 (Dubbs–Edelman). LCB's (κ,λ) pairing
  is reversed in the draft.
- Moment normalizations: D-normalized (Kunisky) vs d-normalized (draft) → the eq. errdef typo.
- GUE ⟨r⟩ benchmark: surmise 0.60266 vs large-matrix 0.5996, used interchangeably.

**Contradictory claims**
- "Gravity/ensemble reproduces every positive moment" vs the draft's own moment table (m₂ 9% off) and its soft tail.
- "Crossing = handle ∝ 1/D²" vs LMRS/Berry (disk crossings O(1), controlled by Δ).
- e^{S0} = D vs LMRS/Berry e^{S0} ∝ BPS count; the draft's own "two entropies" paragraph concedes the mismatch.
- A continuous sub-edge tail vs a bordered-Toeplitz closed form that allows only atoms outside [λ₋, λ₊].
- "Atoms are a resonance near N = 12–13" vs 55.6% at N = 8, 9 and 32.9% at N = 10, 11.
- One-flavor model vs two-flavor monotone tower mixed in the "index" and "chaotic vs deterministic" paragraphs.

**Missing assumptions**
- Haar-typicality of B relative to Π_T (never justified; the random-hyperplane argument needs d ≪ D).
- Absence of harmonic components in the self-energy derivation.
- Choice of R-charge sector / normalization in any gravity comparison.

**Unexplained numerical behavior**
- Why the interior is depleted mid-band by ~2× while "bulk KS" looks small.
- Whether the tail fraction saturates, grows, or is an artifact of the moving threshold λ₋(a(N)).
- The atom excess over rank bounds at N = 13 (+27).
- Even/odd N′ alternation in KS distances.

**Where previous Claude analyses overreached**
- The draft presents a Haar null model as a derived gravity dual and labels conjectures "derived (exact)".
- Machine-precision agreements (Kunisky identity, λ-identity) are presented as verification of physics when they
  are identities.
- Figures and numbers are quoted without code in the repo.
- Citation metadata errors (Chang–Chen–Sia–Yang mis-attributed; missing authors).
- The prior audit is mostly sound but contains its own labelling error (d(N) "interchanged"). It also references a
  "companion research note" with an LMRS second-moment calculation that is not in the repo.

## 8. Recommended project organization (proposal — nothing moved yet)

```
papers/                      primary PDFs (unchanged)
literature/
  notes/                     one note per paper (done)
  synthesis.md, bibliography.bib (done)
research/                    HISTORY (read-only by convention)
  pdfs/, tex/                current draft + prior audit (keep as-is, mark as superseded when rewritten)
  notes/                     dated audit notes (this audit's claim table lives here)
docs/
  project_archaeology.md     this report
  research_questions.md      living
  todo.md                    living
  derivations.md             ONLY verified derivations (block BPS eqs, identities, atom characterization,
                             rank bounds, Haar null model moments/W1) with status labels
  conventions.md             (new) N vs N', a,b,d,D,m, normalizations, benchmarks
src/
  syk_forms.py               (new) single exterior-algebra/BPS-basis module replacing 5 copies
  decoder.py                 O_T, K_T, state-dependent measures, atoms
  nullmodels.py              Wachter/Jacobi (Haar), Weingarten, free Jacobi process
  (existing scripts kept; refactored to import the modules)
scripts/
  run_*.py                   CLI experiment drivers (parameters, seeds, output paths recorded)
  audit/                     verification scripts (done)
tests/                       cross-implementation and exact-identity tests
results/
  data/<date>_<experiment>/  immutable outputs + metadata.json (params, seeds, commit)
  figures/<date>_<experiment>/
paper/                       (new) the rewritten manuscript, separate from research/ history
```

Principles: `research/` is history, `docs/` is current understanding, `results/` is immutable and versioned, and
every figure in a manuscript has a generating script in `scripts/`.

---

## State of the project (starting point for future research)

**What we have.** A well-defined, canonical LMRS-type operator O_T = P_B(1−n_{N′})P_B on the BPS space of one-flavor
N=2 SYK, tied to uplift N → N+1 (Chryssanthacopoulos–Vegh 2026). Its interior level statistics are GUE. Its
spectral density is *close to but measurably not* the Wachter law of two free projectors: mid-band weight moves to
both walls, ~10–15% of states lie below the Wachter edge at N = 12–15, and δm₂ ≈ 0.03–0.04 does not shrink with N.
Exact λ = 0, 1 atoms are dimension-counting effects that vanish for N ≥ 14. For the **Haar null model** we have exact
Weingarten moments and a verified genus-one resolvent W₁.

**What does not hold up.** The current draft's gravitational derivation: the EOW-brane slot, e^{S0} = D, the
trivialized Schwarzian, and "crossing = handle". Also the self-energy formula for O_T, the near-BPS-image explanation
of the tail, the bordered-Toeplitz closed form, and several numerical labels. The Kunisky/Chebyshev "error formula"
is an exact identity (bookkeeping), not evidence.

**Most promising directions.**
1. **LMRS calculation (Q1).** Treat n_{N′} − ⟨n⟩ as a Δ = 1/3 neutral bilinear and compute its zero-energy
   two-point function in super-JT. Test whether the SYK decoder variance decays like N^{−2/3} rather than like the
   Wachter a ~ (√3/2)^N. This is the gravitational calculation the project actually needs.
2. **Reference-subspace structure (Q2).** The enlarged BPS space is a moderate rotation of the LES reference
   B⁰ = B^P_N ⊕ B^{P−1}_N∧e, and the walls/tails are its un-rotated part. Develop perturbation theory in the new
   couplings and a partial-liberation (free Jacobi process) model.
3. **Hygiene first.** Fix the environment, put the repo under git, consolidate the five copies of the exterior-algebra
   code, regenerate the missing data/figures with recorded seeds, and label the Haar model as a null model everywhere.
