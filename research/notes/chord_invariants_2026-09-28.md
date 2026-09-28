# Chord invariants of the uplift decoder (action item 1)

*2026-09-28. Numerical study; all numbers exact linear algebra on disorder realizations, no stochastic trace
estimation. Code: `src/chord_invariants.py` (tested in `tests/test_chord_invariants.py`), driver
`scripts/run_chord_invariants.py`, exact run list `scripts/chord_invariants_campaign.sh`, analysis
`scripts/analyze_chord_invariants.py`. Data: `results/data/chord_invariants_2026-09-28/` (append-only JSONL, one
line per realization; `meta_*.json` record arguments and git commit). Figures:
`results/figures/chord_invariants_2026-09-28/`.*

**RESULTS: §3 (N′ ≤ 15 / 16) and §4 (N′ = 16 / 17 reruns).** Campaign complete 2026-09-28 22:33. The two largest
runs were OOM-killed once and rerun with memory-bounded code under `scripts/memwatch.py`. Theory comparisons:
`docs/derivations.md` D1 (super-Schwarzian) and D2 (finite-λ super-chords, which match r to 2 %).

## 1. Question

A gravitational (chord-diagram) computation of the decoder spectrum would output moments built from
(i) a two-point weight, (ii) crossing weights between distinct operators, and (iii) higher connected pieces. Before
attempting that computation (action item 2), we measure these invariants in SYK, together with their N- and
q-dependence. The three candidate pictures make different predictions:

| picture | single-operator law | crossing weight X | triple crossing (word ijkijk) |
|---|---|---|---|
| commuting / UV (no mixing, perfect uplift) | Bernoulli(b) | 1 | 1 |
| free / Haar P_B (the old draft's model) | Wachter(a, b) | a | X² (free compression) |
| chord / q-Gaussian (DSSYK-type rule) | set by chord weights | X ∈ (0,1) | X³ |

## 2. Definitions and method

Enlarged theory with N′ modes, charge P = ⌊N′/2⌋, one-flavor q-body supercharge, i.i.d. unit complex Gaussian
couplings (same generator and seeds as `src/decoder_fast.py`, `src/q_scan.py`). B = BPS space, d = dim B,
D = C(N′,P), a = d/D, τ = Tr_B/d.

Operators: single-mode slot projectors Π_i = 1 − n_i (Π_{N′−1} is the uplift decoder O_T); k-mode slot projectors
Π_T = ∏_{i∈T}(1 − n_i). With i.i.d. couplings all modes are exchangeable, so the decoder of the "new" mode is
statistically identical to any Π_i. Checked: r(new mode) − r(other modes) is within 2σ at every N′ (q=3, N′=5…14).

Invariants of X = P_B Π P_B (per realization):
- free cumulants κ_n (n = 2…6) of its eigenvalue law, and ρ_n = κ_n / κ_n^{Bern(b)} with b = Tr_D Π / D.
  Free compression (Haar P_B) gives ρ_n = a^{n−1} exactly in the large-D limit (Nica–Speicher). We call
  a_n ≡ ρ_n^{1/(n−1)} the *effective parameters*. "Wachter with renormalized parameter t" means a_2 = a_4 = a_6 = t.
- r ≡ ρ_2 = Var(λ)/[b(1−b)] (fraction of UV variance surviving BPS projection). Here Var uses the realization mean.
- For two distinct operators with Â = P_B(Π − μ)P_B, μ = τ(P Π P): c = τ(ÂᵢÂⱼ)/√(vᵢvⱼ) (covariance),
  F = τ(ÂᵢÂᵢÂⱼÂⱼ)/(vᵢvⱼ) (non-crossing factorization), and **X = τ(ÂᵢÂⱼÂᵢÂⱼ)/(vᵢvⱼ) (crossing weight)**.
  vᵢ = τ(Âᵢ²).
- 6-letter pair words w (letters i, j, k, each twice): w/∏v compared with X^{cr(w)}, cr = number of chord crossings.
  Only the fully crossing word ijkijk (cr = 3) distinguishes free compression (X²) from the chord rule (X³).
- Rigidity fractions η = (SYK − Haar)/(UV − Haar): 0 = free/Haar, 1 = commuting UV. The UV value is computed on the
  full charge sector (r_UV = 1; X_UV from the commuting occupations).

Two exact representations of P_B are used: an orthonormal BPS basis (q=3, small d), and the orthogonal complement
(q=5, where the complement rank r ≪ d). A test checks that they agree to 1e−9 on the same realization; another
checks that the decoder spectrum reproduces `decoder_fast.py`.

**Null model.** A Haar-random projector of the same rank (same D, d), treated identically. It reproduces the
free-compression predictions ρ₄ = a³, ρ₆ = a⁵ and X ≈ a to three digits (e.g. N′=14: ρ₄ 0.077 vs a³ 0.077;
X 0.431 vs a 0.425), which validates the pipeline. Haar KS distances to its own Wachter law are 0.002–0.01.

Sum rule: Σᵢ nᵢ = P on the sector implies c = −1/(N′−1) exactly (observed: −0.25, −0.2, …, −0.071 for
N′ = 5…15). It is not physics, and it contaminates X at O(c²) ≲ 1% for N′ ≥ 10.


## 3. Results (q=3: N′ = 5…15; q=5: N′ = 9…16)

Numbers below are from `scripts/analyze_chord_invariants.py` (rerun 2026-09-28 20:0x on all completed data; the
earlier `summary.csv`/`fits.json` from 19:05 lacked q=3 N′=14–15 and q=5 N′=14–16). Realizations per point: q=3
N′=5–8: 60 SYK / 30 Haar; 9–10: 40/20; 11: 24/12; 12: 16/8; 13: 10/5; 14: 4/2; 15: 3/2. q=5 N′=9–14: 20/10;
15: 8/4; 16: 4/2. Errors are s.e. over realizations. Figures `fig1`–`fig6` in `results/figures/chord_invariants_2026-09-28/`.

Status labels: **[num]** numerical observation at these sizes; **[interp]** interpretation, not established.

### 3.1 The null model validates the pipeline [num]
Haar P_B reproduces free compression to 3 digits for N′ ≥ 7: ρ₄ = a³, ρ₆ = a⁵ (e.g. q=3 N′=15: 0.039 vs 0.039,
0.005 vs 0.005); X_Haar = a + 0.006; F = 1.00; words with cr = 0, 1, 2 equal 1, X, X²; and the fully crossing word
scales as X² (e₃ = 1.94 ± 0.01, constant in N′). The 0.06 offset of e₃ from 2 is a finite-D effect of the null
model itself, and it sets the resolution of the e₃ test.

### 3.2 SYK is not free compression, and the gap grows with N′ [num]
q=3: r (= a₂) exceeds a by an amount that grows with N′: r/a = 1.04 (N′=7), 1.16 (9), 1.27 (11), 1.45 (13),
1.70 (15; r = 0.577 ± 0.003, a = 0.340). A fit over N′ = 7…15 gives r ≈ 0.89 a^0.41 and X ≈ 0.86 a^0.16. Exponential
and power-law fits in N′ are indistinguishable (RSS 1.1e−2 vs 1.2e−2), so the asymptotic form of r (decay to 0
vs a nonzero limit) is **not** determined. Even N′ is systematically more rigid than odd N′ at equal a
(N′=13/14: r = 0.617/0.642). This confirms and extends the audit table in `docs/research_questions.md` Q1.
**Update (same day, `docs/derivations.md` D1):** the LMRS zero-energy two-point function in our single R-charge
sector (j = 0 even N′, −1/6 odd N′) predicts r with no free parameter. It overshoots by 15–24 %, shrinking with N′.
It fits after a *fitted* factor (1 − 2.20/N′), but a plain N′^{−0.44} fits equally well. The parity alternation is
explained, parameter-free, by the j-dependence.

q=5 sits at a = 0.91–0.99 (near-UV). The same effect is present but small: r − a = 0.021 at N′=16.

### 3.3 Single-operator law: Wachter with a renormalized parameter, at the level of moments only [num]
Effective parameters a_n = ρ_n^{1/(n−1)} nearly coincide (q=3 N′=15: a₂, a₄, a₆ = 0.577, 0.590, 0.588; N′=14:
0.642, 0.662, 0.671). So the first six free cumulants look like Wachter(t) with t ≈ r instead of t = a. The first six
moments of Wachter(r) match to 2–6 % (Wachter(a) misses by up to 67 %). There is a systematic upward drift
a₂ < a₄ ≲ a₆ of 3–5 % at even N′, so it is not exact.
**But the distribution is not Wachter(r):** KS distance to Wachter(r) is 0.15–0.23 at q=3 N′ ≥ 11, including
the interior (atoms removed), whereas Haar gives 0.001–0.003. The renormalized-Wachter description is a statement
about low moments/cumulants only.

### 3.4 Pair and word structure: multiplicative up to two crossings, anomalous at three [num]
- Non-crossing factorization holds: F = 0.994–0.995 (q=3), 1.000 (q=5). This is consistent with both the free and
  the chord pictures.
- The crossing weight is O(1) and far above Haar: X = 0.73 (SYK) vs 0.35 (Haar) at q=3 N′=15, decreasing slowly
  with N′ (0.85 at N′=8). Commuting UV: X_UV ≈ 1.
- 6-letter words with cr = 0, 1, 2 obey 1, X, X² to ≲ 1 % (w₂/X² = 0.996–1.005 for N′ ≥ 9).
- **The fully crossing word ijkijk obeys neither rule.** Per-triple e₃ = ln w₃ / ln X (mean ± s.e.):

  | q=3 N′ | 9 | 10 | 11 | 12 | 13 | 14 | 15 |
  |---|---|---|---|---|---|---|---|
  | SYK e₃ | 2.27(1) | 2.29(1) | 2.32(1) | 2.35(1) | 2.35(1) | 2.41(1) | 2.37(2) |
  | Haar e₃ | 1.94 | 1.95 | 1.94 | 1.94 | 1.94 | 1.95 | 1.94 |

  Equivalently, w₃/X³ = 1.17–1.22 (chord rule: 1) and w₃/X² = 0.89–0.94 (free: 1; Haar: 1.03–1.07). e₃ creeps
  upward (odd N′: +0.046, +0.036, +0.016 per step; the increments shrink). Whether it converges to 3, or to a
  value strictly between 2 and 3, cannot be decided at these sizes.
- q=5: e₃ is ill-conditioned there (X ≈ 0.95, so ln X ≈ −0.05). Use ratios instead: w₃/X² = 0.987, w₃/X³ = 1.04
  at N′ = 15–16 (Haar: 1.000, 1.095). SYK sits about a quarter of the way from free toward the chord rule. This is
  qualitatively the same as q=3.

**[interp]** The simplest chord (q-Gaussian) rule, in which each crossing costs an independent factor X, fails
for these operators at accessible N′, and so does free compression. Crossings are not independent: a triple
crossing is less suppressed than free probability predicts and more than the product of three pairwise weights. A
gravitational chord computation would need either an operator- or N′-dependent triple-intersection weight, or
contributions beyond single chord exchange. This is the sharpest discriminating output of the campaign.

### 3.5 Rigidity fractions [num]
η_r (0 = free, 1 = commuting UV), q=3, single mode: 0.03 (N′=5) → 0.24 (9) → 0.30 (11) → 0.33 (13) → 0.36 (15);
even N′: 0.29 (10), 0.34 (12), 0.38 (14). η_X ≈ 0.55–0.59 for N′ ≥ 13. Growth is decelerating; saturation at
η_r ≈ 0.4 is plausible but not established. Heavier probes are less rigid: at N′=15, η_r = 0.36 / 0.29 / 0.22 for
k = 1 / 2 / 3 modes. Rigidity is not a function of a alone. At matched a ≈ 0.9: q=3 N′=6 gives η_r = 0.06, while
q=5 N′=16 gives 0.24 (`fig6`).

### 3.6 What this means for action item 2 (gravity computation)
1. A gravity computation must reproduce **X = O(1)** (0.73–0.85 at q=3), not X ≈ a. This rules out the old draft's
   Haar/Wachter "gravity dual" as a description of SYK at these sizes (consistent with the audit verdict).
2. It must reproduce **F = 1** and multiplicativity up to two crossings. These are satisfied by chord-type rules.
3. It must reproduce **e₃ ≈ 2.3–2.4** (q=3) and its drift. A single-parameter q-Gaussian rule predicts 3; free
   probability predicts 2. This is the main target.
4. The single-operator law should come out as Wachter(t ≈ r) only at the level of the low cumulants.

## 4. Memory incident and reruns (2026-09-28)

The last lines of groups C (q=3 N′=16, `--method basis`) and E2b (q=5 N′=17) were OOM-killed while running
concurrently on the 32 GB machine; they left only `meta_q3_N16_*.json` / `meta_q5_N17_*.json` and no raw rows.
Causes: (i) a full `np.linalg.svd` of the dense 8736×12870 constraint matrix, with estimated peak ≈ 16 GB; (ii) an
unbounded cache of r×r matrices (r ≈ 3060, 150 MB each) keyed on raw bytes, with about 100 entries.
Fixes in `src/chord_invariants.py`: null space of H = M†M = {Q, Q†} via partial `eigh`, with a certified spectral
gap (used only for D > 7000, so no earlier run is affected); a byte-bounded LRU cache keyed on letter multisets;
and `estimate_memory_gb` / `--dry-run`. Regression check: new code on stored realizations (q=3 N′=12 seed 0 SYK
and Haar, q=5 N′=13 seed 0) reproduces all 495 recorded numbers to ≤ 4e−16. The `eigh` path agrees with `svd` to
4e−14 (tests `test_eigh_nullspace_matches_svd`, `test_tiny_cache_gives_same_words`). Reruns: campaign groups F1,
F2, under `scripts/memwatch.py` (caps 13 + 11 GB within the 30 GB budget; peaks logged in `memwatch.jsonl`).

**Rerun results (complete 22:33; peaks 7.54 GB / 6.98 GB against caps 13 / 11 GB; `memwatch.jsonl`).**
q=3 N′=16: 2 SYK + 1 Haar realizations, 4 modes, 1 triple. q=5 N′=17: 2 SYK + 1 Haar, 4 modes, 1 triple.
Summaries and figures regenerated.

| q | N′ | a | r | X | a₂, a₄, a₆ | e₃ (per triple) | η_r (k=1,2,3) | η_X (k=1) | Haar r, X |
|---|---|---|---|---|---|---|---|---|---|
| 3 | 14 | 0.425 | 0.642 ± 0.004 | 0.769 ± 0.003 | 0.642, 0.662, 0.671 | 2.408 ± 0.009 | 0.377, 0.293, 0.214 | 0.594 | 0.425, 0.432 |
| 3 | 16 | 0.340 | 0.601 ± 0.005 | 0.755 ± 0.002 | 0.601, 0.626, 0.636 | 2.42 ± 0.03 | 0.396, 0.316, 0.244 | 0.625 | 0.340, 0.346 |
| 5 | 16 | 0.913 | 0.934 | 0.949 | 0.934, 0.936, 0.937 | w₃/X² = 0.987 | 0.242, 0.162, 0.104 | 0.411 | 0.913, 0.914 |
| 5 | 17 | 0.874 | 0.908 | 0.931 | 0.908, 0.909, 0.908 | w₃/X² = 0.979, w₃/X³ = 1.052 | 0.267, 0.193, 0.148 | 0.451 | 0.874, 0.874 |

What the new points change:
- Every trend in §3 continues. None reverses.
- **η_r, even N′:** 0.340 → 0.377 → 0.396 (N′ = 12, 14, 16). The increments are shrinking (0.037, 0.019), which
  supports but does not prove saturation near 0.4–0.45. q=5 is still rising (0.214 → 0.242 → 0.267).
- **e₃ (q=3):** 2.42 ± 0.03 at N′=16. Its upward creep is consistent with §3.4, but the error bar at N′=16 is too
  large to resolve it. It is still clearly between 2 and 3.
- **r:** agrees with the finite-λ super-chord prediction (0.589; D2) to 2 %, and with the pre-registered 1/N′-fit
  range 0.597–0.604.
- **Haar null:** again reproduces free compression exactly (r = a, X ≈ a + 0.006, e₃ = 1.94).
- Statistics at N′ = 16/17 are thin (2 SYK realizations). Treat them as trend confirmation only.
