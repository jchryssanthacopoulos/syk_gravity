# Resonant chords, layer 1: the first-order dressing spectrum (2026-10-26)

Short record. The derivation and the results table are in `docs/derivations.md` D11. Context: D7 (central chords) and
D10 §4 (the window-edge resonance). The plan for the resonant chord problem has three layers:
1. first order (this note);
2. the finite-N hard edge;
3. the self-dressing R (continuum → continuum chords).

## Code
- `src/chord_dressing.py` (backward compatible; `tests/test_chord_dressing.py` passes unchanged):
  - `DressingFamily(om, fam, inner='svd')`: inner energies from singular values. D7's Gram cut drops the near-zero edge
    modes at N = 15.
  - `old_model_from_uplift(OldModel(F))`: admissible sources in any sector.
  - `chord_population(..., F=)`, `chord_t1(..., F=)`, `analyze_family(..., Fsec=)`. The "chord_rescaled" layer is
    skipped when s_in ≤ ½, since t₁ is then infinite.
- `scripts/run_chord_resonance.py`: resonant h family at ℓ′ = −2, plus the central control ℓ′ = 0 (both families).
- `scripts/analyze_chord_resonance.py`: error budget, tails, `fig_layer1`.
- Test: `tests/test_uplift_window.py::test_dressing_family_in_resonant_sector`.

## Runs (under `scripts/memwatch.py`; outputs in `results/data/chord_resonance_2026-10-26/`)
| job | command | seeds | etas | peak | time |
|---|---|---|---|---|---|
| small | `run_chord_resonance.py --Np 8 10 12 14` | 9000–9002 | 24 | 0.11 GB | 142 s |
| N16 | `run_chord_resonance.py --Np 16 --no-control --neta 8` | 9000, 9001 | 8 | 3.35 GB | 1232 s |

Couplings are as in `results/data/uplift_window_2026-10-26`, so the own-η spectra equal the window runs' σ² (checked to
2e−15 at N′ = 12).

## Results (D11 §3)
- The central control reproduces D7's error budget, which validates the pipeline.
- Resonant sector, chord vs exact first order:
  - W +8.5, +6.7, +4.8 % (N′ = 12, 14, 16), as good as central;
  - T₂ −8.8, −7.5, −7.0 %, about twice central, mostly from the chord continuum;
  - tail amplitude P·√K(100) low by 26, 21, 12 %. That deficit is entirely continuum: the planar model with the
    exact S matches the exact tail to ≤ 2.5 %.
  - Every error shrinks with N′.
- The Kronecker model with the measured BPS heterogeneity fails in resonant sectors from N′ = 12 (T₂ 0.71 vs 0.52). The
  heterogeneity is dominated by a few edge modes.

## Layer 3: the back-reaction operator R and the exact resonant decoder (same day; D11 §4)
Code:
- `scripts/run_chord_layer3.py`: models for R inside the exact self-dressed h-block pencil:
  - exact; Haar-rotated; flat spectrum; c_R/E_w;
  - full chord (chord clouds + chord R);
  - diag R (D6 mean field); R = 0 (first order).
- `scripts/analyze_chord_layer3.py`: table and `fig_layer3`.

Runs (under memwatch; outputs in `results/data/chord_layer3_2026-10-26/`, see the naming note below):
- small: `--Np 10 12 14 --seeds 9000 9001 9002 --haar 3`, peak 0.25 GB, 75 s.
- N16: `--Np 16 --seeds 9000 9001 --haar 2`. Seed 9000 finished (peak 4.84 GB, 833 s); seed 9001 was interrupted by
  a session restart before writing, and rerun alone (`log_N16_seed9001.txt`, peak 4.94 GB, 821 s).
  The first log is kept as `log_N16_seed9000_then_interrupted.txt`.

Findings:
- R = A A† with A = P_Co E H₀^{−1/2} on coexact(ℓ′−5). It has no resonant denominators, and its rank is
  dim coexact(ℓ′−5): r_R/n = 0.113, 0.172, 0.236 at N′ = 10, 12, 14, against exp(−10λ/9) = 0.108, 0.162, 0.215.
- The soft edge modes have 74–90 % of their weight in ker R. R is far from diagonal in the edge basis (diagonal fraction
  0.11–0.17), so the diagonal mean field fails.
- Free (Haar) model: R → URU† reproduces the exact h-block decoder (T₂ to ≤ 0.5 %, P(λ < 0.1) to 2–5 %, N′ = 10–16).
  So do the flat spectrum (rank + mean only) and the c_R/E_w spectrum. The deepest tail P(λ < 0.03) is underestimated,
  increasingly with N′ (0.012 vs 0.021 at N′ = 16): the exact range R is slightly aligned with the softest edge modes.
- Full chord model:
  - T₂ is off by −1.4, −7.9, −7.7, −8.8 % at N′ = 10, 12, 14, 16, inherited from layer 1's chord clouds;
  - the R-induced shift of T₂ is reproduced (−0.009 / −0.014 / −0.021 against exact −0.005 / −0.018 / −0.022);
  - the floor is ~10⁻², against 10⁻⁶ at first order.

Naming note: the layer-3 directories were first created as `*_2026-10-27` by mistake (the runs are dated 2026-10-26).
They were renamed after the runs finished, and the scripts and D11 updated. The meta files inside keep the original
absolute paths.
