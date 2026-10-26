# The uplift anywhere in the window (2026-10-26)

Short record. The derivation is `docs/derivations.md` D10, and the question is Q13 in `docs/research_questions.md`.
The report is `research/tex/progress_2026-10-26.tex`. The literature input is `literature/notes/kanazawa_wettig_2017.md`.

## What was done
- **Literature.** Kanazawa–Wettig 2017 was read (§6, Apps. A–B) and noted.
  - App. B (N^±_f, N^z_f) is reproduced by our rank computations for N_c = 5, 9, 12, 13 (one realization each).
  - Dictionary: KW's V⁺ is coexact, KW's V⁻ is exact. CV App. A uses the opposite labels.
- **Code** (backward compatible: `tests/test_uplift.py`, `tests/test_chord_dressing.py` unchanged and passing).
  - `src/uplift.py`:
    - `OldModel(..., F=)` works in any sector. It computes the connecting maps δ_h, δ_g and the admissible sources
      VF, VF1; the full harmonic bases are kept as BF, BF1.
    - `harmonic_basis` handles any degree. It uses the same draws as `bps_basis_fast`, and the output was checked
      identical.
    - Also added: `uplift_bookkeeping`, `strip_walk_count`, `bps_dim_q3`, `wedge_pairing`, `connecting_form`,
      `OldModel.dyson`, and safe `pinv` / `pinv_apply` for empty maps.
    - `decoder_blocks` now handles an empty family. It used to crash on the overlap of an empty matrix.
  - `src/chord_dressing.py`: `chord_sectors` and `chord_jump_law` take an optional F.
  - Scripts:
    - `scripts/run_uplift_window.py` (tasks `sectors`, `forms`);
    - `scripts/analyze_uplift_window.py` (checks, tables, figures, chord tail, Hill estimates);
    - `scripts/plot_uplift_window_map.py` (jump maps).
  - `tests/test_uplift_window.py`: 7 tests.
- **Couplings.** rng = default_rng(seed), c = random_couplings(rng, C(N′, q)). All sectors of one (N′, seed) share c.
  At q = 3 and half filling these are the realizations of `results/data/uplift_2026-10-12` (seeds 9000, 9001).

## Runs (all under `scripts/memwatch.py`; outputs in `results/data/uplift_window_2026-10-26/`)
| job | command (`scripts/run_uplift_window.py …`) | seeds | peak | time |
|---|---|---|---|---|
| A | `--task sectors --p 3 --Np 5 … 14` | 9000–9002 | 0.83 GB | 919 s |
| B | `--task forms --p 3 --N 4 … 14`; `--task forms --p 5 --N 6 … 14` | 9000, 9001 | 0.50 GB | 38 s |
| C | `--task sectors --p 5 --Np 6 … 14` | 9000, 9001 | 1.51 GB | 1048 s |
| D | `--task sectors --p 3 --Np 15 --F 6 7`; `--Np 16 --F 7 8` (ℓ′ ≤ 0; ℓ′ > 0 follows by mirror symmetry) | 9000, 9001 | see log | see log |

Incidents (no effect on the results):
- The first launches of A and C (logs `log_sectors_p3_small_failed1.txt`, `log_sectors_p5_failed1.txt`) crashed in
  `decoder_blocks` on their first one-family sector, before writing any record. After the fix both were relaunched.
- A and D started in the same second. D's meta file `meta_sectors_p3_20261026T103527.json` overwrote A's failed-run
  meta. Meta file names now include the PID.
- Two smoke-test meta files, whose records went to the scratchpad, were moved out of the results directory.

## Results
See D10. In brief:
- Every check passes:
  - dimensions, connecting-map ranks against the predictions, uplift residuals ≤ 6.4e−13, closed-form vs direct
    decoder ≤ 2.6e−11, and mirror symmetry ≤ 2.7e−12;
  - 139 records at q = 3 and 114 at q = 5.
- Exceptional states are the radical of the connecting form: radicals 1, 3, 1 at N = 4, 8, 12 (q = 3) and 1, 1, 0 at
  N = 6, 10, 14 (q = 5).
- q = 5 at N′ = 19: Q is injective on Λ²(19) (rank 171, σ_min/σ_max = 0.89, two realizations). By CV's parity, an odd
  number of exceptional states follows, beyond KW's tested range.
- Resonant sectors |ℓ′| = q − 1:
  - edge-multiplet gap 6e−6 to 8e−7 at N = 11–13;
  - first-order dressing tail α = 0.43 / 0.48 (N′ = 12 / 14) against the chord value ½; amplitude 0.66 vs 0.50 at
    K = 100;
  - the exact decoder shows no wall pile-up.
- T₄ = T₂ − H_first holds to ≤ 0.009 in all extensive sectors at N′ = 9–14.

## Scratch (not in the repo)
`scratchpad/stage3/`:
- `sweep.py`, `resonance_peek.py`, `dyson_check*.py`, `form_check.py`, `radical_eta.py`: exploratory versions. Their
  numbers agree with the production runs where both exist.
- `radical_eta.py`: the N = 8 radical has dimension 3 for every η tried. Its subspace moves with η (principal-angle
  cosines 0.05–0.49).
