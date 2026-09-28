# Audit scripts (2026-09-28)

Independent checks written during the project archaeology audit (`docs/project_archaeology.md`).
They verify or falsify specific claims in `research/tex/gravitational_dual_fragility.tex`.
None of them writes to `results/`.

**Environment.** Run in a clean venv (tested: Python 3.12, numpy 2.5.3, scipy 1.18.1, sympy 1.12). The
Anaconda base environment on this machine is broken for scipy/matplotlib (NumPy 1.x/2.x ABI mismatch).

| Script | Claim tested | Result (2026-09-28) |
|---|---|---|
| `weingarten_W1_check.py [kmax]` | Haar genus expansion: disk moments, m_k^{(1)} (draft eqs. m2handle–m4handle, m₅^{(1)}), W₁ closed form (eq. W1) | Exact Weingarten via class-algebra inversion. k = 1…5: no odd 1/D terms; disk = Wachter; genus-1 = draft formulas = W₁ series. **All True.** At a=b=½: m₂ = 3/8 − ε²/8 (draft's "−¼D⁻²" is wrong). Runtime ~minutes for k=5. |
| `kernel_identity_check.py` | Draft Sec. 7: Tr O^k = D m_k^W + Σ K_{k,ℓ} T_ℓ | Holds to ~1e−14 with **d**·m_k^W for Haar, random-coordinate and strongly aligned (non-free) projector pairs, k ≤ 7 ⇒ an algebraic identity, not a physics test. With the literal D·m_k^W it fails by O(D). |
| `reference_subspace_check.py` | Exploratory (research Q2): principal angles between B^P_{N′} and B⁰ = B^P_N ⊕ B^{P−1}_N∧e; weight of edge eigenvectors in B⁰ | N′ = 10, 11, 12; 2 seeds each (seeds 0, 1 of `decoder_fast.couplings`). Mean cos² = 0.81–0.88; eigenvectors with λ<0.05 have 97–99.6% weight in B^{P−1}_N∧e; λ>0.95 ~99% in B^P_N; bulk (0.2<λ<0.8) ~51–58%. **Preliminary** (small N, few seeds; includes exact atoms). |

Smoke test of existing code (same venv, outputs in a scratch directory, not committed):
`src/decoder_fast.py --N 8 9 --seeds 2` and `src/protected_atoms.py --N 8 9 10 11 --seeds 1` run correctly.
Atom counts (λ=1, λ=0, fraction of d): N=8: 26, 19, 55.6%; N=9: 45, 45, 55.6%; N=10: 44, 36, 32.9%;
N=11: 80, 80, 32.9%. Per-realization m₁ at N′=9 is 0.536–0.560 (b = 0.556); at N′=10, m₁ = 0.5 exactly.
