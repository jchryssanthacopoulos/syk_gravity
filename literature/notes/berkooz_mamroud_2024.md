# Berkooz & Mamroud (2024) — A cordial introduction to double scaled SYK

**Citation.** M. Berkooz, O. Mamroud, *A Cordial Introduction to Double Scaled SYK*, arXiv:2407.09396 (v2, Nov 2024).
File: `papers/berkooz_mamroud_2024.pdf`. (Read: Secs. 1–2, 4 (rules), 6.2, 7.2; Secs. 3, 5 skimmed.)

## Main content
- Double-scaling limit p, N → ∞ with λ ∝ p²/N fixed (eq. 1.2). Moments of H become sums over **chord diagrams**;
  each crossing of two H-chords weighs q = e^{−2p²/N} (Sec. 2.1). Transfer matrix / q-deformed oscillator (Sec. 2.2).
- **Matter** (Sec. 2.3): probe operators M_Δ of length p̃ with independent Gaussian couplings (eq. 2.26), Δ = p̃/p.
  An H-chord crossing an M-chord weighs **q̃ = e^{−2pp̃/N} = q^Δ** (eq. 2.28). Two-point function eq. (2.30).
  Crossed four-point function via an R-matrix / 6j symbol of U_q(su(1,1)) (eqs. 4.1–4.5). The general rule is
  that crossing of chords of lengths p₁, p₂ weighs e^{−2p₁p₂/N}.
- Triple-scaling limit → Schwarzian / JT (Sec. 3.1); quantum-group structure and non-commutative AdS₂ (Sec. 5).
- **Supersymmetric DSSYK** (Sec. 6.2): N=1 via chords between Q insertions; N=2 with oriented chords reducible to
  the ordinary transfer matrix. Used by Boruch–Lin–Yan (JHEP 12 (2023) 151) to relate the fraction of ground
  states to chord-Hilbert-space wavefunctions.
- Sec. 7.2: beyond strict double scaling, sample-to-sample 1/N^p effects (global spectral fluctuations) are much
  larger than the exponentially small ramp/plateau contributions. Their gravity interpretation is open.

## Relationship to our project
- Supplies the concrete **finite-Δ crossing rule** missing from the draft. The draft asserts crossings are
  suppressed by 1/D² (the Haar rule). In chord language they are suppressed by q^{Δ₁Δ₂}-type weights depending on
  operator size. For our UV operator n_{N′} (size p̃ = 2, fixed), the crossing weight with other fixed-size chords
  e^{−2·2·2/N} → 1 at large N: the **light-operator regime where crossed and uncrossed diagrams have equal weight**
  (matching LMRS eq. 169 and Berry p. 60). Freeness (non-crossing dominance) is not expected for Π_T.
- Suggests a computational route for research Q1/Q4: supersymmetric chord techniques (Boruch–Lin–Yan) for the
  BPS-projected two-point function of a fixed-size bilinear at finite N.
- The 1/N^p sample-to-sample effects are a possible source of finite-N non-universality in decoder moments, and
  should be distinguished from genuinely non-Haar structure.
