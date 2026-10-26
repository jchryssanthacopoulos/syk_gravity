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
