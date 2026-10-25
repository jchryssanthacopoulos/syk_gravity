# Chords for the two-flavor uplift: black holes dress, gravitons do not (2026-10-25)

Short record. The full account is `research/pdfs/progress_2026-10-25.pdf` (source `research/tex/progress_2026-10-25.tex`)
and `docs/derivations.md` D9. Question: Q12, Stage 2, in `docs/research_questions.md`.

## What was done
- **Literature check.** The papers in `papers/` contain no chord treatment of the two-flavor model: CCSY §3.3 and HTZ §3
  use G–Σ only. We derived the two-flavor chord rules by redoing BBNR §§4.1–4.5 with two flavors (D9 §1).
- **Code.**
  - `src/two_flavor_chords.py`: the chord dictionary (λ, s, δ, D(s), continuum edges, two-point functions, uplift
    bookkeeping).
  - `scripts/run_two_flavor_chords.py`: tasks window, blocks, resonance, pairing, tower, edgegap, twisted.
  - `scripts/analyze_two_flavor_chords.py`: summary and figures.
  - `tests/test_two_flavor_chords.py`: 6 tests, run as `.venv/bin/python tests/test_two_flavor_chords.py`.
- **Couplings.** As in `scripts/run_two_flavor_uplift.py`: C1 = random_{symmetric,generic}_couplings(N+1,
  default_rng(seed)), old model = restriction. So old model N, seed s is the same realization as in
  `results/data/two_flavor_2026-10-23`. The window task draws its N'-pair model at size N' with the same rule.

## Runs (all under `scripts/memwatch.py`; outputs in `results/data/two_flavor_chords_2026-10-25/`)
| job | command (`scripts/run_two_flavor_chords.py …`) | seeds, models | peak | time |
|---|---|---|---|---|
| small | `--tasks blocks tower resonance pairing --N 3 4 5 6`; `--tasks window --Np 4 5 6 7` | 1–3, both | 0.19 GB | 392 s |
| N7 sym | `--tasks tower blocks resonance pairing --N 7` | 1–2, sym | 2.55 GB | 2687 s |
| N7 generic | `--tasks blocks resonance pairing --N 7` | 1–2, generic | 2.70 GB | 2636 s |
| window | `--tasks window --Np 8`; `--Np 9 10 --dmax 3000` | 1–2, both | 2.11 GB | 1463 s |
| edgegap | `--tasks edgegap --N 3 5 7` (1–3); `--N 9` (1–2) | both | 6.80 GB (generic N = 9 in parallel: 6.76, 6.78 GB) | 3062 s (+1187, 1172 s) |
| twisted | `--tasks twisted --N 4 5 6 7` | 1–3, both | 0.25 GB | 100 s |

A first window attempt (N' = 8, old eigensolver) was stopped by hand after 538 s (peak 1.75 GB) before writing any record
(its meta file is `meta_sym_window_20261025T164429.json`).

The N = 9 edge-gap diagonalizations take about 30 min each on one core, so the generic N = 9 seeds were run in parallel:
- the main edgegap job (`log_edgegap.txt`) ran N = 3, 5, 7 for both models and N = 9 symmetrized seeds 1–2;
- its generic N = 9 step (`meta_generic_edgegap_20261025T174347.json`) was stopped seconds after it started, before writing
  anything;
- generic seeds 1 and 2 ran as separate jobs (`log_edgegap_gen9s1.txt`, `log_edgegap_gen9s2.txt`), appending to the same
  `edgegap_N9_generic.jsonl`.
The dry-run estimates (4.9 GB at N = 7) were conservative.

Code changes during the campaign (no change to results):
- `Spectra.get` was optimized after the small-N window runs. The nonzero spectrum now comes from the Gram matrices of
  the two Q maps instead of a third diagonalization of the Laplacian.
  - Checked against the earlier version at N' = 7, seed 1: identical BPS counts; m1 and gaps agree to 1.7e−12.
  - The first N' = 8 window attempt, with the old code, was stopped before writing any record and rerun.
- The edgegap and twisted tasks were added while the N = 7 resonance runs were in progress. Those runs use code
  without them, and the other tasks were not changed.

## Results
See D9 §5 and the report. In brief:
- **Dictionary.** Two-flavor chords = one-flavor N=2 chords at λ = 10/N, s = (2δN_ψ − δN_χ)/5. E, K, Z are chiral with
  δ = 3/5, 4/5, 2/5. Twisted chords (symmetrization) are 1/N on black holes and act as a V insertion.
- **Window.** Every BPS sector at N' ≤ 8 lies in |s| ≤ ½ except the outer tower sectors. At N' = 9 only tower states lie
  outside. (1,1) at N' = 9, 10 has the tower (sym) and nothing (generic).
- **Dictionary accuracy.** Obstruction weights 1.05–1.21 (chord/exact); continuum integrals 0.80–1.08; continuum edges
  0.65–0.97, all with log-scatter ≤ 0.35 (edges ≤ 0.47).
- **Resonance.** At N' = 8, J-resonant bulk sectors with a large edge target have mean κ 125 / 48 and largest
  7.8×10⁴ / 1.5×10⁴ (sym / generic), against 1.0 and 30 non-resonant. Same-charge sectors with small boundary targets
  show no enhancement (a control).
- **Edge gap.** Falls fast to N = 7 and stalls at N = 9 (3.3–4.8×10⁻⁵ of the mean, both models, seeds 1–2).
- **Tower.** Currents cancel to ≤ 2.5e−16. Its individual emissions are black-hole-sized, also at the resonance.
  E–K coherence of black holes is the same in both models; the tower's cos² = 1.
- **Pairing.** Fortuitous c within [−0.023, 0.021]; tower c = n(N'−n)/N'².

## Scratch (not in the repo)
- `scratchpad/stage2/blocks*.py`, `angles.py`, `gaps.py`, `tower_resonance.py`, `slotcorr.py`: exploratory versions of
  the tasks.
- Their numbers agree with the production runs where both exist.
- `angles.py` (principal angles between the E- and K-emission subspaces, without removing the tower) turned out to be
  dominated by dimension counting at N ≤ 6 (cos θ = 1 forced). It is not used. The `twisted` task replaces it.
