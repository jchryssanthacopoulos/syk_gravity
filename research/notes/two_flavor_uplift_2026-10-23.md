# The two-flavor pair uplift: monotone vs fortuitous (2026-10-23)

Short record; the full account is `research/pdfs/progress_2026-10-23.pdf` (source `research/tex/progress_2026-10-23.tex`)
and `docs/derivations.md` D8. Question: Q12 in `docs/research_questions.md`.

## What was done
- Model: symmetrized two-flavor N=2 SYK (CCSY eq. 3.18; CV §6.2). C = (G + G^{T_jk})/2, G i.i.d. complex Gaussian
  (`random_symmetric_couplings`). Control: generic couplings (`random_generic_couplings`).
- Uplift N → N+1 adds the pair (ψ_e, χ_f) at the top two bits; the old model is the restriction of the same couplings.
  Slots T = (t_ψ, t_χ).
- Code: `src/two_flavor.py` (model), `src/two_flavor_uplift.py` (`PairUplift`), `tests/test_two_flavor.py` (16 tests,
  run as `.venv/bin/python tests/test_two_flavor.py`).

## Runs (all under `scripts/memwatch.py`)
| run | command | seeds | peak |
|---|---|---|---|
| symmetrized N = 2–6 | `run_two_flavor_uplift.py --N 2 3 4 5 6 --seeds 1 2 3 4 5 --model sym` (+ `--tasks e2`) | 1–5 | 0.44 GB |
| generic N = 2–6 | `run_two_flavor_uplift.py --N 2 3 4 5 6 --seeds 1 2 3 --model generic --tasks dims cv bare first` (+ `e2`) | 1–3 | 0.52 GB |
| symmetrized N = 7 | `run_two_flavor_uplift.py --N 7 --seeds 1 --model sym` (dims, cv, bare completed; stopped during `first`, see below) | 1 | 4.59 GB (cap 5) |
| symmetrized N = 7 | `run_two_flavor_uplift.py --N 7 --seeds 1 --model sym --tasks first tower` (rerun with per-sector cache release) | 1 | 3.20 GB (cap 8) |
| symmetrized N = 7 | `run_two_flavor_uplift.py --N 7 --seeds 1 --model sym --tasks e2` | 1 | 0.23 GB |
| symmetrized N = 7 | `run_two_flavor_uplift.py --N 7 --seeds 2 --model sym --tasks cv bare tower e2` | 2 | 3.65 GB |
| generic N = 7 | `run_two_flavor_uplift.py --N 7 --seeds 1 --model generic --tasks cv bare e2` | 1 | 3.27 GB |

The first seed-1 N = 7 job ran the code of commit de2a4d9. That code kept every sector's first-order pseudo-inverses in
memory, and the job reached 4.59 GB of its 5 GB cap during `first`. Records are written only when a task completes, so
nothing was written for `first`, and no data were lost. It was stopped by hand and `first tower` was rerun with code that releases those caches after
each sector.

Code versions: the `meta_*.json` files record `de2a4d9` (the N ≤ 6 runs and the first seed-1 N = 7 job) or
`de2a4d9-dirty` (the `e2` runs, seed 2 and generic at N = 7, and the seed-1 rerun).
- The dirty states differ from de2a4d9 by the added E₂ method and `e2` task, and, for the seed-1 rerun only, by the
  cache release. Neither changes what the other tasks compute.
- Both were committed unchanged in d69dbb2 (`git diff d69dbb2 -- src scripts/run_two_flavor_uplift.py` is empty).
- Later commits change only the analysis and plotting script.

Outputs: `results/data/two_flavor_2026-10-23/{dims,cv,bare,first,tower,e2}_N*_{sym,generic}.jsonl`, `meta_*.json`
(arguments, git state), `log_*.txt`, `memwatch_*.jsonl`. Analysis: `scripts/analyze_two_flavor_uplift.py` →
`summary.json` and `results/figures/two_flavor_2026-10-23/fig_{cv_fidelity,tower_costs,bare_counts,first_order}.{pdf,png}`.

Exploratory scratch scripts (not in the repo) also used seeds 300+N, 500+N, 600+N, 700+N, 900+N, 950+N, 4000+N. Their
numbers agree with the production runs. One exception: seed 705 (N = 5) has an unusually small old-model gap (0.0067),
which inflates first-order dressing costs by ~10×. It is not used anywhere in the report.

## Results (details and status tags in D8)
- **Exact.**
  - Block form Q_{N+1} = Q₀ + ψ_e†E + Kχ_f + ψ_e†Zχ_f.
  - [Q₀,E] = [Q₀,K] = 0 and [E,K] = {Q₀,Z}.
  - Protection identity [E,V] = K (symmetrized only).
  - Eight slot equations.
  - Theorem A (one mode with obstruction), and the two-step dimension formula.
  - The tower is bare, with slot weights (N+1−n)/(N+1) and n/(N+1).
  - The first-order formulas.
- **Numerical laws** (0 exceptions outside N = 5):
  - four-slot injectivity;
  - rk D_T = min(d_new, R_T);
  - dim(B ∩ im D_T) = max(0, b + rk D_T − R_T);
  - κ_T(tower) = CV bound ⟺ ker D_T = 0;
  - bare = generic count + tower;
  - BPS dimensions are identical for symmetrized and generic couplings.
- **Corrections to working assumptions made earlier in this session.**
  - The pair spectral sequence does *not* degenerate at E₂: d₂ ≠ 0 in four sectors per N.
  - "Monotone = coherent, fortuitous = incoherent" is false: fortuitous states interfere too, and the tower is
    distinguished by exact cancellation.
- **Unexplained.** Three one-dimensional generic-position failures at N = 5 (near-central sectors, both models, all
  seeds). They are not tolerance artifacts: the old-model gaps are ≥ 0.03.

## Revision 2026-10-25: closed forms (D8 §9; report §9)
- **Runs** (code commit 813a905, under `scripts/memwatch.py`):
  - `run_two_flavor_uplift.py --tasks closed` for symmetrized seeds 1–5 (N = 2–6), generic seeds 1–3 (N = 2–6), and seed 1
    of each model at N = 7. Output: `closed_N*_{sym,generic}.jsonl`.
  - `--tasks levels` (commit f586878), same seeds at N = 2–6. Output: `levels_N*.jsonl`.
  - Peaks: 0.43 GB (N ≤ 6), 5.50 GB (N = 7).
- **Pure-dressing counterexamples**: `scripts/find_pure_dressing_states.py`
  - `--model generic --N 4 --seed 2 --sector 3 3`
  - `--model sym --N 4 --seed 1 --sector 2 2`
  - Output: `pure_dressing_*.json`.
- **Exact results.**
  - Cohomology = H(old BPS data, d₁ + d₂), with d₂ = P_H(Z − EQ₀⁺K + KQ₀⁺E)P_H (homological perturbation).
  - BPS states = (1 + 𝒦)⁻¹h with explicit admissibility.
  - Integrating out the middle slots gives Dyson equations whose self-energies are EGE† + KGK† and K†GK + E†GE, with a
    leak through the middle cohomology.
  - Uniqueness criterion ρ < 2.
  - Slot weights and decoder pencils.
- **Numerical results.**
  - All closed-form checks pass to ≤ 7e−13, with zero dimension mismatches.
  - Four-slot injectivity fails on a codimension-2 set of couplings (pure-dressing BPS states, both models). Earlier
    statements that it "holds" are superseded: it is generic, not universal.
- Scratch checks (not in the repo): `scratchpad/closed/cf1.py`–`cf5.py`, `gscan.py`, `g2scan*.py`.
