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
  - `tests/test_uplift_window.py`: 9 tests (incl. KW App. B rows for N = 5, 9, 13 and the resonant chord weight).
- **Couplings.** rng = default_rng(seed), c = random_couplings(rng, C(N′, q)). All sectors of one (N′, seed) share c.
  At q = 3 and half filling these are the realizations of `results/data/uplift_2026-10-12` (seeds 9000, 9001).

## Runs (all under `scripts/memwatch.py`; outputs in `results/data/uplift_window_2026-10-26/`)
| job | command (`scripts/run_uplift_window.py …`) | seeds | peak | time |
|---|---|---|---|---|
| A | `--task sectors --p 3 --Np 5 … 14` | 9000–9002 | 0.83 GB | 919 s |
| B | `--task forms --p 3 --N 4 … 14`; `--task forms --p 5 --N 6 … 14` | 9000, 9001 | 0.50 GB | 38 s |
| C | `--task sectors --p 5 --Np 6 … 14` | 9000, 9001 | 1.51 GB | 1048 s |
| D | `--task sectors --p 3 --Np 15 --F 6 7`; `--Np 16 --F 7 8` (ℓ′ ≤ 0; ℓ′ > 0 follows by mirror symmetry) | 9000, 9001 | see below | see below |
| E | `--task sectors --p 3 --Np 16 --F 7 --seeds 9001 --out sectors_p3_N16_seed9001.jsonl` (parallel to D) | 9001 | see below | see below |

Job D was **stopped by hand** (rc 143, SIGTERM) right after its seed-9000 N′ = 16 F = 8 record. Its next item,
seed 9001 F = 7, duplicated job E, and the remaining F = 8 seed 9001 duplicates D6. Peak 7.82 GB, 3653 s. N′ = 16 records:
F = 7 (ℓ′ = −2) seeds 9000 (job D, 1531 s) and 9001 (job E); F = 8 (ℓ′ = 0) seed 9000 (job D, 1500 s; T₂ = 0.6075 =
D6's exact value).

**A numerical failure caused by the resonance itself (N′ = 16, seed 9001, ℓ′ = −2).**
- The first job-E record had d = 2187 + 1 instead of 2187 + 0, with split overlap 0.99999975. The dimension check
  flagged it.
- Cause: `harmonic_basis` decided ranks by Gram eigenvalues σ² > 1e−10 σ²_max (as `bps_basis_fast` does).
  - At N = 15 the smallest mode of the window-edge multiplet (Λ⁶ → Λ⁹) fell below that threshold for this seed,
    so a coexact edge mode was counted as harmonic at level −3.
  - `pinv`, which thresholds σ itself at 1e−10 σ_max, kept the same mode.
- Fix: `harmonic_basis` now takes ranges from SVDs with the `pinv` criterion. All earlier records passed the
  dimension checks, so their harmonic spaces are unaffected (the bases agree up to round-off).
- The bad record is kept as `superseded_gramrank_N16_seed9001.jsonl` (log `superseded_gramrank_log_N16_seed9001.txt`),
  outside the analysis glob. It was rerun as `sectors_p3_N16_seed9001_svdrank.jsonl`.
- Diagnosis (`scratchpad/stage3/diag_n16.py`, peak 0.81 GB): on Λ⁶(15) → Λ⁹(15), the Gram eigenvalues below
  1e−10 · max number 454 for seed 9000, which is exactly the exact 6-forms, C(15,3) − 1. For seed 9001 they number
  455: the extra one is the near-zero edge mode.

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
  - 142 records at q = 3 (N′ ≤ 16) and 114 at q = 5. The superseded N′ = 16 seed-9001 record is excluded; its
    rerun (`sectors_p3_N16_seed9001_svdrank.jsonl`, peak 5.42 GB, 1261 s) passes.
- Exceptional states are the radical of the connecting form: radicals 1, 3, 1 at N = 4, 8, 12 (q = 3) and 1, 1, 0 at
  N = 6, 10, 14 (q = 5).
- q = 5 at N′ = 19: Q is injective on Λ²(19) (rank 171, σ_min/σ_max = 0.89, two realizations). By CV's parity, an odd
  number of exceptional states follows, beyond KW's tested range.
- Resonant sectors |ℓ′| = q − 1:
  - edge-multiplet gap (seed means) 5.9e−6, 1.6e−6 at N = 11, 13; 8.0e−8 and 1.5e−10 at N = 15 (two seeds);
  - first-order dressing tail, Hill α (K > 100) = 0.43 ± 0.06, 0.48 ± 0.04, 0.48 ± 0.03 at N′ = 12, 14, 16, against
    the chord value ½;
  - amplitude P·√K at K = 100: 0.66, 0.66, 0.59 vs chord 0.50, 0.51, 0.52;
  - the exact decoder shows no wall pile-up.
- T₄ = T₂ − H_first holds to ≤ 0.009 in all extensive sectors at N′ = 9–14.

## Scratch (not in the repo)
`scratchpad/stage3/`:
- `sweep.py`, `resonance_peek.py`, `dyson_check*.py`, `form_check.py`, `radical_eta.py`: exploratory versions. Their
  numbers agree with the production runs where both exist.
- `radical_eta.py`: the N = 8 radical has dimension 3 for every η tried. Its subspace moves with η (principal-angle
  cosines 0.05–0.49).

## Paused 2026-10-26 (at the user's request: battery)
- The rerun of N′ = 16, seed 9001, F = 7 (`sectors_p3_N16_seed9001_svdrank.jsonl`) was stopped by hand before writing a
  record.
- To resume:
  1. Relaunch it (about 30 min, peak about 6 GB, memwatch cap 10 GB).
  2. Run `scripts/analyze_uplift_window.py`. The tail statistics now use the h family at ℓ′ = −2 only.
  3. Replace the placeholders in `research/tex/progress_2026-10-26.tex`: N16ROWS, N16ALPHA, N16PEAK, COMMIT.
  4. Update the N′ = 14 rows of the D10 §5 table and of report Table 1 to the final 3-seed values: ℓ′ = −2 is
     0.477 → 0.455 → 0.455, H 0.133/0.142, error −0.0092; ℓ′ = 0 is 0.566 → 0.619 → 0.636, error −0.0034.
  5. Update the means (N′ = 14 first-order mean 851, exact 787; gap 1.6e−6) and the Hill n (mirror double count
     removed).
  6. Build the PDF and commit.
- Resumed the same day. The rerun completed: d = 2187 + 0 as expected, T₂ = 0.4303, peak 5.42 GB, 1261 s. Its aborted
  first attempt is logged in `log_sectors_p3_N16_seed9001_svdrank_aborted_pause.txt`.
