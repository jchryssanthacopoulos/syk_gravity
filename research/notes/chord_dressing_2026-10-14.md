# The first-order dressing spectrum from chords (2026-10-14)

Short record; the full account is `research/pdfs/progress_2026-10-14.pdf` (source `research/tex/progress_2026-10-14.tex`)
and `docs/derivations.md` D7.

- **Closed forms** (`src/dssyk_n2.py`):
  - BPS→continuum two-point function of a charged chiral single-chord operator = normalized Al-Salam–Chihara measure
    (`bps_continuum_weight`); point mass = BLY 4.14 (`bps_bps_weight`); time integral = ₂φ₁
    (`time_integrated_two_point`).
  - Verified against the transfer-matrix engine (`src/chord_transfer.py`) in `tests/test_chord_dressing.py`.
- **BBNR sign:** the eigenvector coefficients of (5.10) belong to +cos φ in (5.3). Use H_n(−cos φ|q²) when combining with
  (5.3).
- **Finite-N tests** (`scripts/run_chord_dressing.py`, `results/data/chord_dressing_2026-10-14/decomp_N*.jsonl`):
  - chord W^(1)/exact = 1.04–1.09 (N′ = 10–16), no fit;
  - planar-law error +3–5 % (soft lower edge), chord inputs +1–3 %, factorization < 1 %;
  - crossed diagrams ≤ 0.5 % of E tr 𝒟³ (`crossing_N*.jsonl`).
- **Heterogeneity** (`scripts/analyze_bps_heterogeneity.py`, `heterogeneity_N*.jsonl`):
  - sum rule B_c + Σ_i P E_i†GE_iP = F P;
  - corr(𝒟̄, B_c) = −0.73 … −0.93;
  - coherent share of the IR variance 0 → 0.66 (N′ = 10–15).
- **Pre-registered, not run:** `prereg_N17.json` (committed in 097cd6e before any N′ = 17 computation).
- **Stopped when the report was finalized:** the second N′ = 16 realization of the decomposition, and the N′ = 16
  heterogeneity analysis.
