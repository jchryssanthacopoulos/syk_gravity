# Haikin, Zamir & Gavish (2017) — Random subsets of structured deterministic frames have MANOVA spectra

**Citation.** M. Haikin, R. Zamir, M. Gavish, *Random Subsets of Structured Deterministic Frames have MANOVA
Spectra*, arXiv:1701.01211 [cs.IT]. File: `papers/haikin_zamir_gavish_2017.pdf`. (Read: abstract and Sec. 1.)

## Main content
- Empirical study: for many deterministic equiangular / tight frames, a random k-subset of frame vectors has
  Gram spectrum indistinguishable from Wachter's MANOVA law (and the Jacobi ensemble at finite size).
- Proposed universal fluctuation law |E_K Ψ(λ) − Ψ(f^MANOVA)|² = C n^{−b} log^{−a} n with universal exponents (eq. 1).
- Purely numerical; conjecture later proven (weak convergence) by Kunisky (2023).

## Relationship to our project
- Precedent that *structured, non-Haar* projection pairs can nevertheless have MANOVA spectra. So Wachter
  agreement by itself does not imply the BPS space is Haar-typical, and Wachter disagreement is not automatic
  for structured subspaces.
- Their fluctuation-scaling methodology (convergence rate vs n) is a model for testing whether SYK decoder
  deviations shrink with N (they apparently do not, over the accessible range).
