# The uplift map, the decoder walls, and what localization is

*2026-10-12. Follows `research/notes/decoder_T4_fourarc_2026-10-05.md` (D5: decoder T₄ = liberation + isotropic +
localization). Derivations: `docs/derivations.md` D6. Code: `src/uplift.py` (tests `tests/test_uplift.py`),
`scripts/run_uplift_walls.py`, `scripts/run_uplift_blocks.py`, `scripts/analyze_uplift_heterogeneity.py`,
`scripts/plot_uplift_mechanism.py`. Data: `results/data/uplift_2026-10-12/`; N′ = 15 OU/decoder data in
`results/data/ou_fourarc_2026-10-05/N15_{a,b}/`. Figures: `results/figures/uplift_2026-10-12/`. Polished write-up: `research/pdfs/progress_2026-10-13.pdf`.*

Labels: **[derived]** proven here or in D6; **[num]** numerical observation; **[interp]** interpretation;
**[pre-registered]** prediction written down (with timestamp) before the test was run.

## 0. Summary

**Question.** D5 split the decoder's non-free fourth moment into free liberation, a small isotropic part, and a
"localization" excess (decoder − isotropic OU at equal T₂; 0.03–0.047 in T₄ at N′ = 10–14). What is localization,
and can it be predicted?

**Answer.**
1. **The N′-mode BPS space is an explicit, exact deformation of the (N′−1)-mode one [derived; verified to 10⁻¹⁵].**
   - Ψ(h, g) = (h + x, g + z) for old harmonic states h ∈ B_o^F and g ∈ B_o^{F−1}, where (1 + 𝕃)(x, z) = (εK′g, −εKh).
   - 𝕃 is *anti-Hermitian*. It is built from the old Green's operator acting on the new mode's pair-creation operator
     E = η∧.
   - Consequences:
     - an exact orthogonal splitting B = B_h ⊕ B_g;
     - a closed-form 2×2 block decoder, whose diagonal blocks are functionals of the old model's integrated pair
       two-point functions 𝒟 = E†H₀⁻¹E and 𝒟′ = EH₀⁻¹E†;
     - monotone self-dressing.
2. **Walls are arithmetic [derived + pre-registered].** Exact wall states are Lefschetz kernels (η∧h = 0). The
   near-wall hard edge exists only while dim(exact forms of degree F+2) ≲ dim B_o^F. It was predicted from old-model
   data to collapse at N′ = 15, and it does.
3. **Localization is not the walls [pre-registered prediction refuted].** L(N′=15) = 0.0454 ± 0.0003, flat across
   N′ = 13–15, although the walls collapse.
4. **Localization is a heterogeneity deficit, and the old model predicts it [identity + num].**
   - T₄ = T₂ − H with H = ⟨y(1−y)⟩, so L = H^OU − H^dec exactly.
   - The first-order uplift spectrum of the old model reproduces H^dec to 1–7 % (N′ = 11–15). So T₄ is predicted to ≤ 0.008
     (liberation alone: 0.03–0.06), and L to 2–25 %.
   - Out-of-sample N′ = 16 [pre-registered]: pre-registered H residual window [0, +0.005]; measured +0.0030 and +0.0042 (2 seeds), so T₄ is predicted to +0.003 / +0.004. Pass. Side predictions on wall fractions and on the y < ½ fraction narrowly missed (D6 §6).

## 1. Exact structure (summary of D6)

- **Split.** Λ^F(N′) = X ⊕ Y∧e with X = Λ^F(o) and Y = Λ^{F−1}(o). Then ω = ω₀ + η∧e, Q₀ = ω₀∧ and E = η∧, with
  [E, Q₀] = 0 (η is a 2-form).
- **BPS equations.** Q₀a = 0, Q₀†b = 0, Q₀b = −εEa and Q₀†a = εE†b, with ε = (−1)^F.
- **Hodge pieces.** a = h + x (x ∈ Ex_X), b = g + z (z ∈ Co_Y). K = Q₀⁺E, K′ = (Q₀†)⁺E†, K_H = K|_{H_X},
  K′_H = K′|_{H_Y}, L = K|_{Ex_X}.
- **Adjoint relation [derived].** P_Ex K′ P_Co = L†, from [E, Q₀] = 0 and [E†, Q₀†] = 0. Hence
  (1 + 𝕃)(x, z) = (εK′_Hg, −εK_Hh) with 𝕃 = [[0, −εL†], [εL, 0]] anti-Hermitian, always invertible.
- **Orthogonal splitting [derived].** The uplifts of B_o^F (b coexact) and of B_o^{F−1} (a exact) are orthogonal,
  because (1+𝕃)⁻¹†(1+𝕃)⁻¹ = diag((1+L†L)⁻¹, (1+LL†)⁻¹) is block-diagonal. Measured overlaps ≤ 1.1·10⁻¹⁵,
  N′ = 9–15.
- **Canonical uplift [derived; added 2026-10-13].** The BPS uplifts of an old state h are Ψ(h, g) for any g, and the freedom
  in g is CV's coset ambiguity. By the splitting, Ψ(h, 0) is the unique minimum-norm uplift.
- **Closed-form decoder [derived; reproduces the N′-mode decoder to 10⁻¹⁰].** 𝒩v = μ𝒢v with
  - 𝒢 = diag(1 + K_H†(1+R)⁻¹K_H, 1 + K′_H†(1+R′)⁻¹K′_H);
  - 𝒩_hh = 1 − K_H†(1−R)(1+R)⁻²K_H;
  - 𝒩_gg = −1 + K′_H†(1−R′)(1+R′)⁻²K′_H;
  - 𝒩_hg = −2εK_H†L(1+R′)⁻²K′_H;
  - R = LL†, R′ = L†L.
- **Dressing cost [derived].** σ² ∈ spec P E†H₀⁻¹E P = spec ∫₀^∞dτ P E†e^{−τH₀}E P.
  - The diagonal blocks are functionals (by push-through) of 𝒟 = E†H₀⁻¹E on the closed F-forms and
    𝒟′ = EH₀⁻¹E† on the coclosed (F−1)-forms.
  - The h–g block also needs K′|_{H_Y} itself.
- **Monotonicity [derived].** Ordered self-dressed h-eigenvalues are ≥ the first-order 1/(1+σ²), and g-eigenvalues
  are ≤ σ′²/(1+σ′²). Mixing adds 2‖C₂‖²/d ≥ 0 to T₂.

## 2. Walls [derived; num]

- **Exact atoms.** λ = 1 iff η∧h = 0, and λ = 0 iff ι_η g = 0.
- **Refined bound.** Per wall, b − (exact-form dimension). Totals:

  | N′ | 10 | 11 | 12 | 13 | 14 | 15 |
  |---|---|---|---|---|---|---|
  | measured (identical across realizations) | 90 | 80 | 160 | 41 | 82 | 0 |
  | refined bound | 90 | 79 | 158 | 27 | 54 | 0 |

  The excess is unexplained.
- **Hard edge.** Ratio r = e_{F+2}/b_F: 0.44, 0.67, 0.99, 0.96, 1.31, 1.33 at N′ = 10, 12–16.
- **[pre-registered, confirmed]** The fraction of eigenvalues with min(λ, 1−λ) < 0.001:
  - first order from old-model data: 0.117 (N′ = 14) → 0.003–0.012 (N′ = 15);
  - exact: 0.128 → 0.012–0.023 (6 realizations).
- **Near-wall states are first-order exact.** Upper quantiles agree to ≤ 3·10⁻³. Wall eigenvectors have ‖h‖² ≈ 0.99
  and ‖η∧h‖/‖h‖ ≈ 1.1 (bulk: 4–8).

Wall fractions, ε = 0.1 / 0.03 / 0.01 / 0.001 (`walls_N*.jsonl`):

| N′ | exact | first order |
|---|---|---|
| 13 | 0.55 / 0.37 / 0.26 / 0.123 | 0.49 / 0.32 / 0.23 / 0.11 |
| 14 | 0.58 / 0.39 / 0.27 / 0.13 (3 seeds) | 0.50 / 0.34 / 0.24 / 0.117 (13 seeds) |
| 15 | 0.49–0.54 / 0.28–0.31 / 0.15–0.18 / 0.012–0.023 | 0.42–0.47 / 0.22–0.27 / 0.11–0.14 / 0.003–0.012 |
| 16 | 0.541–0.551 / 0.318–0.322 / 0.180 / 0.021 (2 seeds) | 0.461–0.465 / 0.264–0.267 / 0.144–0.146 / 0.010–0.011 (same seeds) |

## 3. Localization at N′ = 15: the wall hypothesis fails [pre-registered; num]

- **Prediction (written before the run).** If localization were the wall pile-up, L would drop at N′ = 15 with the
  hard edge.
- **Run.**
  - `run_ou_fourarc.py --Np 15 --K 8 --samples 4 --seed0 {0,100} --fast-basis --pair-starts 2 --parity-sets 1`, run as
    two memwatch jobs with a 6 GB cap each (peak 3.5 GB).
  - K = 8 instead of 16 for time (about 15 min per sample). The lag-interpolation error of K = 8 was checked at
    N′ = 12–14: ≤ 2·10⁻⁴ in T₄^OU(x*).
- **Result (decoder s = 1).**
  - T₂^dec = 0.5754, T₄^dec = 0.4524, x* = 0.721 ± 0.001 (8 samples).
  - T₄^OU(x*) = 0.4070, T₄^lib = 0.3918.
  - **L = 0.0454 ± 0.0003**, against 0.0432 (N′ = 13) and 0.0467 (N′ = 14).
  - Shares of T₄^dec − T₄^free: liberation 77 %, isotropic 5.8 %, localization 17 % (N′ = 14: 78 / 3 / 19 %). Output: `results/data/ou_fourarc_2026-10-05/analysis_output_N15.txt`.

## 4. The heterogeneity prediction [num]

**Method.**
- For each realization, first order from the old model alone: y₁ = ((1−σ²)/(1+σ²))² (h-states) and the analogue for
  σ′ (g-states). From these, H_first = ⟨y₁(1−y₁)⟩.
- 13 seeds at N′ = 8–14, 6 at N′ = 15, 2 at N′ = 16 (plus the 3 first-order seeds behind the pre-registration).
  Seeds 9000 + r; couplings from `random_couplings`, the same law as all project data.
- Compare with H^dec = T₂^dec − T₄^dec (all-site means of the OU campaign), and on the same realizations with the
  exact decoder.

| N′ | T₄^dec | T₂ − H_first | error | liberation error | OU(x*) error | L measured | L predicted | H_first | H^dec |
|---|---|---|---|---|---|---|---|---|---|
| 8 | 0.7778 | 0.7884 | +0.0105 | -0.0151 | -0.0175 | 0.0175 | 0.0280 | 0.0291 | 0.0396 |
| 9 | 0.6666 | 0.6694 | +0.0029 | -0.0312 | -0.0293 | 0.0293 | 0.0322 | 0.0625 | 0.0653 |
| 10 | 0.6790 | 0.6869 | +0.0079 | -0.0315 | -0.0321 | 0.0321 | 0.0400 | 0.0587 | 0.0666 |
| 11 | 0.5799 | 0.5793 | -0.0006 | -0.0444 | -0.0380 | 0.0380 | 0.0374 | 0.0897 | 0.0891 |
| 12 | 0.5964 | 0.6023 | +0.0059 | -0.0446 | -0.0417 | 0.0417 | 0.0475 | 0.0838 | 0.0897 |
| 13 | 0.5101 | 0.5077 | -0.0025 | -0.0539 | -0.0432 | 0.0432 | 0.0408 | 0.1107 | 0.1083 |
| 14 | 0.5283 | 0.5323 | +0.0040 | -0.0541 | -0.0467 | 0.0467 | 0.0507 | 0.1045 | 0.1085 |
| 15 | 0.4524 | 0.4490 | -0.0034 | -0.0606 | -0.0454 | 0.0454 | 0.0420 | 0.1264 | 0.1230 |

(OU-campaign all-site decoder means; H_first is the mean over 13 seeds, 6 at N′ = 15; `heterogeneity_analysis.json`.)

**Notes.**
- The per-seed difference H_exact − H_first is nearly deterministic (SE ≤ 0.0008). It alternates with parity:
  +0.0079, −0.0012, +0.0045, −0.0026, +0.0034, −0.0041 for N′ = 10–15.
  - At even N′ it shrinks.
  - At odd N′ it grows in magnitude. This is the main caveat for extrapolation.
- **Parameter-free variant.** T₂ from the BLY chord formula (D2) together with H_first: errors in T₄ are +0.025, +0.013, +0.018, +0.003, +0.010, −0.006, +0.002, −0.012 at N′ = 8–15 (≤ 0.012, i.e. ≤ 2.7 % of T₄, for N′ ≥ 11). These errors are dominated by the 1–2 % error of the chord T₂ (D2). No input from the N′-mode model is used.

## 5. Mechanism [num; derived pieces in D6 §5]

- **First order → self-dressed → exact (block decomposition).**
  - T₂ rises by 0.026–0.052 (self-dressing; monotone by Courant–Fischer) and then by 0.008–0.019 (h–g mixing,
    = 2‖C₂‖²/d).
  - H moves by ≤ 0.003 under self-dressing (N′ ≥ 11), and by +0.004 (even N′) or −0.001 to −0.003 (odd N′) under
    mixing. Full table in D6 §5; data in `blocks.jsonl`.
- **Per-state tracking** (`plot_uplift_mechanism.py`, N′ = 14, 2 seeds, `mechanism_N14.json`).
  - Projected first-order states keep 97.4 % of their norm in the exact BPS space (min 84 %).
  - E[y|j] − y₁ ≈ +0.11 for y₁ < 0.6, falling to +0.006 above 0.9.
  - E[y(1−y)|j] − y₁(1−y₁) runs from +0.070 (y₁ < 0.1) to −0.045 (y₁ ≈ 0.65).
  - The mean cancels to +0.003. Figure: `fig_uplift_mechanism.pdf`.
- **Mean-field dressing.** Replace R by ρ_j = ⟨k_j,Rk_j⟩/‖k_j‖² along each first-order dressing direction. This
  reproduces the self-dressed T₂ to ≤ 0.005 (N′ = 11–14). Reading: the first-order dressing cloud (n_e = 1) of an
  h-state partially converts, with weight R(1+R)⁻¹, into exact components of a (n_e = 0), pushing the state toward the
  wall.
- **Withdrawn** (an earlier reading in this session): "first order misassigns states by y ↔ 1 − y". The folded
  distributions agree (KS 0.02–0.03 at N′ ≥ 13), but tracking shows a shift plus spread, not a reflection.
- **Also falsified this session.** Scalar and amalgamated free unitary Brownian rotations of B⁰ (scratch models)
  fall on the self-liberation curve (N′ = 10: T₄ = 0.6515 vs 0.6519 at T₂ = 0.75). The SYK excess over that curve
  grows like (1 − T₂)² along the site-coupling family and saturates at about +0.031 (T₄), +0.05 (T₆), +0.062 (T₈).

## 6. Interpretation [interp]

- The decoder angle of a BPS state is set by its **dressing cost**. This is the Euclidean-time-integrated
  zero-temperature two-point function of the new mode's pair-creation operator (a random two-fermion matter
  operator) in the old BPS sector. The τ-integral is cut off by the gap of the charge-(F+2) sector.
- At equal T₂ the decoder's angles are more polarized than those of an isotropic deformation:
  - less weight at intermediate y ≈ ½;
  - more weight at the extremes, y ≈ 1 (λ near 0 or 1) and y ≈ 0 (λ ≈ ½).

  In the uplift picture y is a function of the dressing cost (y = ((1−σ²)/(1+σ²))² at first order). So localization
  is a statement about how dressing costs are distributed across old BPS states, and first order captures this
  distribution well enough to give H.
- The gravitational counterpart is the spectrum (not just the trace) of an integrated matter two-point function in
  the BPS sector. A BLY-type matter-chord computation with an H⁻¹ insertion is the concrete next step.

## 7. Status and open items

- Exact and verified:
  - the uplift formula, the adjoint relation, the anti-Hermitian system, the orthogonal splitting, the closed-form
    block decoder, monotone self-dressing (D6 §§1–5);
  - wall atoms = kernels;
  - T₄ = T₂ − H.
- Numerical, robust over N′ = 10–15 (+16 out of sample):
  - H_first ≈ H^dec, hence the T₄ / localization prediction;
  - the arithmetic hard edge.
- Open:
  - why H is preserved (approximate cancellation);
  - the odd-N′ growth of the residual;
  - the atom excess;
  - the chord computation of the 𝒟, 𝒟′ spectra.

## 8. Reproducibility

```
# first order (+ exact with --exact), per realization: walls_N{N'}.jsonl
.venv/bin/python scripts/run_uplift_walls.py --Np 14 --nreal 3 --seed0 9000 --exact
.venv/bin/python scripts/run_uplift_walls.py --Np 14 --nreal 10 --seed0 9003          # first order only
# N' = 15, 16 under the watchdog (peaks: 3.0 GB at N'=15; 10.8 GB at N'=16)
.venv/bin/python scripts/memwatch.py --limit-gb 13 --log results/data/uplift_2026-10-12/memwatch_walls.jsonl -- \
    .venv/bin/python scripts/run_uplift_walls.py --Np 16 --nreal 1 --seed0 9000 --exact
# first order -> self-dressed -> exact decomposition: blocks.jsonl
.venv/bin/python scripts/run_uplift_blocks.py --Np 10 11 12 13 14 15 --nreal 2 --seed0 9000
# prediction table + figure
.venv/bin/python scripts/analyze_uplift_heterogeneity.py --Np 10 11 12 13 14 15 16
# mechanism figure
.venv/bin/python scripts/plot_uplift_mechanism.py --Np 14 --seeds 9000 9001 --json results/data/uplift_2026-10-12/mechanism_N14.json
# N' = 15 OU campaign (two jobs; merge afterwards)
.venv/bin/python scripts/memwatch.py --limit-gb 6 -- .venv/bin/python scripts/run_ou_fourarc.py --Np 15 --K 8 \
    --samples 4 --seed0 0 --fast-basis --pair-starts 2 --parity-sets 1 --out results/data/ou_fourarc_2026-10-05/N15_a
cat results/data/ou_fourarc_2026-10-05/N15_{a,b}/ou_p3_N15_F7.jsonl > results/data/ou_fourarc_2026-10-05/ou_p3_N15_F7.jsonl
```
Pre-registrations: `results/data/uplift_2026-10-12/prereg_N16.json` (with its first-order inputs
`prereg_inputs_walls_N16_first.json`). The N′ = 15 wall-collapse and localization-drop predictions were stated in the
session log before the runs. Their outcomes are in D6 §4 and §6.
