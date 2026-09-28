# Review of `scaling_limits.md`

*Reviewed 2026-09-28. The reviewed note is unchanged. Legend: ✅ verified · ⚠️ correct but needs qualification ·
❌ unsupported or incorrect.*

## Verdict

The **kinematics is correct and useful**: exact b, the ϑ₄ formula for a in the double-scaled limit, η* ≈ 0.685, and
rank-forced atoms for η > η*. The **dynamical verdicts are not established**:
- "double-scaled = pure Wachter" assumes the freeness it needs to prove;
- the Airy↔Bessel labelling of η* is wrong;
- the "triple-scaled" section is inconsistent with the note's own correction box.

## 1. Claims that check out

| Claim | Status | Notes |
|---|---|---|
| b = m/D = 1 − p/(N+1) | ✅ | Exact. b = ½ requires N′ = N+1 even. For odd N′ at P = ⌊N′/2⌋, b = (N′+1)/(2N′) (e.g. 0.538 at N′=13), so "b = ½ identically" is asymptotic. |
| At b = ½: u = ½, λ± = ½ ± √(a(1−a)), p₀ = p₁ = max(0, 1 − 1/(2a)) | ✅ | |
| General p₀ = max(0, 1 − b/a, 1 − (1−b)/(1−a)) (Setup) | ⚠️ | Spurious third term. Per-BPS normalization: p₀ = max(0, 1 − b/a), p₁ = max(0, 1 − (1−b)/a). |
| a = d/D → ϑ₄(0, e^{−2η}), η = q²/N | ✅ (DS limit) | Derivation correct: the two alternating image sums in D̃_p − D̃_{N−p−q} (CV eq. 125) combine into the two-sided theta series. Conditional only on CV Conj. A.1 (rank formula) and the local CLT. Numerical check below. |
| η* = 0.685 from ϑ₄ = ½ | ✅ | |
| Single-scaled: a ~ e^{−r(q)N}, r(q) = −ln cos(π/2q), r(3) = 0.1438; Wachter band → δ(λ−½) | ✅ | |
| Atoms of weight 1 − 1/(2a) survive at large N for η > η* | ✅ / ⚠️ | Forced by rank: d > m ⇒ dim ker ≥ d − m. Calling them "fortuity's protected directions" overstates it: they are dimension-counting atoms, present equally in the Haar model. |

### Numerical check of the ϑ₄ formula

Exact a from CV eq. 125 at the central sector p = N/2 versus ϑ₄, along fixed η (N = enlarged size):

| η | ϑ₄(0, e^{−2η}) | q=3 | q=5 | q=9 | q=15 | q=25 | q=41 |
|---|---|---|---|---|---|---|---|
| 0.3 | 0.0749 | 0.0617 (N=30) | 0.0681 | 0.0734 | 0.0744 | 0.0746 | 0.0748 (N=5604) |
| 0.685 | 0.5001 | 0.4248 (N=14) | 0.4974 | 0.4979 | 0.4998 | 0.5000 | 0.5000 |
| 1.0 | 0.7300 | 0.6429 (N=10) | 0.7001 | 0.7210 | 0.7268 | 0.7288 | 0.7296 |

Convergence to ϑ₄ is clear, but not uniform in q. At q = 3 and the sizes we can diagonalize, ϑ₄ overestimates a by
4–9% (N = 8, 10, 12, 14: exact 0.771, 0.643, 0.526, 0.425 vs ϑ₄ 0.789, 0.671, 0.559, 0.459). The fixed-q
asymptotics of ϑ₄ give rate π²/(8q²) = 0.137 at q = 3, not the exact 0.144. The note's first line ("every a checked
against eq. 125 and ED") should state which formula was used where.

## 2. Claims that do not hold

### 2.1 "Double-scaled ⇒ pure Wachter" (§2) — ❌ assumed, not derived
- The argument is that "overlap chords have crossing weight 1/D² = e^{−2S0}, which vanishes, so the ensemble is at
  the free point". But 1/D² crossing weights **are** the Haar/Weingarten model: they presuppose that P_B is free
  (Haar-rotated) relative to Π_T, which is the claim to be proven. The same issue is flagged in the draft audit
  (`gravitational_dual_fragility_audit.md` §2.6–2.7, §4.3).
- "Provable (both parameters O(1), Weingarten + eq. 125 give a Wachter + handles series)" is true only of the Haar
  model. Eq. 125 fixes the dimension d, not the relative position of the BPS space and the slot.
- What double scaling does establish: *if* the decoder is free, the resulting Wachter law is non-degenerate
  (O(1) band width). That is worth stating; it is not the same claim.

### 2.2 The chord heuristic points the other way — ⚠️ open
In double-scaled SYK, a probe of fixed size p̃ crosses supercharge/Hamiltonian chords with weight
q̃ = q^Δ, Δ = p̃/q (Berkooz–Mamroud eq. 2.28). For Π_T = 1 − n_{N′}, p̃ = 2, so Δ = 2/q → 0 and q̃ → 1. That is a
**light** probe, the regime where crossed and uncrossed diagrams have equal weight (LMRS eq. 169; Berry p. 60).
It is the opposite of the non-crossing/free point. Physically, a fraction ~q/N → 0 of supercharge terms involves
mode N′, so n_{N′} becomes approximately conserved. That tends to push O_T's spectrum toward {0,1} (the "UV
imprint"), not toward a Wachter band. This is a heuristic, not a result. The defensible statement is that in the
DS limit both r ≡ Var(λ)/[b(1−b)] and a are O(1), and **whether r(η) = a(η) is the open question**.

### 2.3 "t₂ → 0 is non-perturbative in 1/D, invisible to the genus series" (§2 caveat) — ❌
t₂ ≈ δm₂ = m₂^SYK − m₂^W is an O(1) shift of the second moment, visible at leading order, not a non-perturbative
effect. At q = 3 it **grows** with N: δm₂ = 0.006 → 0.057 over N′ = 7…14, with r/a rising 1.03 → 1.53
(`docs/research_questions.md` Q1; `scripts/audit/variance_scaling.py`). It is the central open question, not a caveat.

### 2.4 "η* is an Airy↔Bessel edge-universality transition" (§2) — ❌
At b = ½ the lower band edge is λ₋ = ½ − √(a(1−a)). It is strictly positive, with square-root (soft) behaviour, on
**both** sides of a = ½. For a > ½ the atoms sit at λ = 0 separately, below a soft edge. The edge touches 0 and
becomes hard (density ~ λ^{−1/2}, Bessel kernel) only at a = b exactly, i.e. in a window where d − m = O(1)
(cf. Collins 2004 hard-edge regime). η* is a critical point, not a boundary between an Airy phase and a Bessel phase.

### 2.5 "Triple-scaled" limit (§3) — ❌
- **Not DSSYK triple scaling.** In DSSYK, triple scaling zooms the *energy spectrum of H* onto the Schwarzian
  regime. Here the zoom is onto λ → 0 of the decoder spectrum. The link to H goes through E, which the note's own
  correction box says is a Rayleigh quotient (tail E ~ 10⁻² vs smallest nonzero H_N eigenvalue ≈ 2.25; occupied
  component ≳ 99.99% old-BPS harmonic). So "resolves the near-extremal throat", "the individual near-BPS states,
  the gap, ρ_H" and ρ_tail dλ → ρ_H(E)dE as a density of states do not follow. §3 should be rewritten around the
  correction, not alongside it.
- **Change of variables.** ρ_tail(λ) = σ_a²/(1−λ)² ρ_H(σ_a²λ/(1−λ)) is algebraically right *if* σ_a² is constant.
  σ_a² is state-dependent (draft: 55.0 ± 3.0 at N = 12), so this is an approximation.
- **σ_a² ~ C(N/2, q−1) ~ N²/8.** The second step holds only for q = 3. In the DS limit (q ~ √N) C(N/2, q−1) is not
  ~N²/8. Even at q = 3, the measured σ_a² = 55 at N = 12 is ~3× N²/8 = 18. The absolute scale also depends on the
  coupling normalization (unit variance in the code vs FGMS's 2J/N² vs CV's J(q−1)!/N^{q−1}). Only E/σ_a² is
  meaningful.
- **"ρ ~ λ^{−0.6}"** disagrees with repo data: fitted exponents of about −0.23 to −0.11 at N = 14–15
  (`results/data/softtail_scaling.csv`, column rho_exp); the draft reports −0.4 → −0.13 over N = 12–15.
- **Modular-duality reading** ("band = winding/chord rep, tail = momentum/energy dual; triple scaling = passage to the
  dual frame"): the ϑ₄ modular transform relates the small-η and large-η expansions of the *count* a, not of the
  decoder spectrum. Unsupported; remove or mark as speculation.

### 2.6 Minor
- §1: "fortuity's protected core is a vanishing finite-N resonance". The atom fraction decays monotonically
  (55.6% at N=8,9 → 32.9% at 10,11 → 5.6% at 12,13 → 0 at 14–16) and tracks rank bounds; it is not a resonance
  (draft audit item 1.13).
- §1: "the entire nontrivial spectrum is the near-BPS tail" at fixed q. Not what the data show: r ≈ 0.65 at N′ = 14
  (while a = 0.43), i.e. a large O(1) spread persists. What SYK converges to at fixed q is open (Q1).
- Notation: the note mixes N (old) and M = N+1 (enlarged) within the count formula. Fix to the conventions in
  `docs/project_archaeology.md` §7 (N old, N′ = N+1).

## 3. Suggested corrected summary

> b → ½; in the double-scaled limit η = q²/N fixed, a → ϑ₄(0, e^{−2η}) (proven modulo CV Conj. A.1), so the
> Wachter law would be non-degenerate, with rank-forced atoms of weight 1 − 1/(2a) for η > η* = 0.685. At fixed q,
> a → 0 and the Wachter band collapses. **Whether the SYK decoder is actually free in the DS limit, i.e. whether
> r(η) = a(η), is open.** The light-probe chord rule (Δ = 2/q → 0) suggests not; at q = 3 the deviation from
> Wachter grows with N.

## 4. Test that would settle §2

Compare r − a at matched a across q, moving toward the DS limit. For a ≈ 0.77–0.80: q = 3 at N′ ≈ 7–8 (dense ED,
done: r/a ≈ 1.03–1.06) versus q = 5 at N′ ≈ 20–22 (needs the matrix-free estimator `src/q_scan_mf.py`; dense ED
is infeasible there). If the note is right, r − a should shrink as q grows at fixed a. Under the light-probe
heuristic, it should not.
