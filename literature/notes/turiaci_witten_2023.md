# Turiaci & Witten (2023) — N=2 JT supergravity and matrix models

**Citation.** G. J. Turiaci, E. Witten, *N = 2 JT Supergravity and Matrix Models*, arXiv:2305.19438 (v4, Nov 2023).
File: `papers/turiaci_witten_2023.pdf`. (Read: Secs. 1, 2.1–2.3, 2.5; Secs. 3–7 and appendices skimmed — 126 pp.)

> **Citation note.** The current draft cites this as "G. J. Turiaci, E. Witten *et al.*". There are two authors.

## Main question
Which random-matrix ensemble is holographically dual to N=2 JT supergravity on surfaces of arbitrary topology?

## Main results
1. **Ensemble** (Sec. 2.3, Fig. 2): supermultiplets of different R-charge are **statistically independent**. Each
   (k, k+q̂) multiplet sector is an Altland–Zirnbauer (α,β) = (1+2ν, 2) ensemble for the supercharge restricted
   to that sector, with ν the number of BPS states. Time-reversal variants in Sec. 6.
2. **BPS counts do not fluctuate**. The numbers L_k⁰ of BPS states are fixed by the disk and receive no wormhole
   fluctuations (Sec. 2.3, Sec. 3.2). BPS wavefunctions are random (Sec. 1).
3. Leading spectrum from the disk (Sec. 3.1): e^{S0} BPS states in a restricted range of R-charges, a gap for most
   (q̂, δ), then a continuum (Fig. 1 summarizes N=0,1,2,4).
4. Gravity side: measure on the moduli space of N=2 hyperbolic surfaces (eq. 4.52), an N=2 Mirzakhani-type
   recursion for volumes (eq. 5.31), and a proof that the volumes satisfy the matrix-model loop equations (Sec. 5.4;
   App. B handles the logarithmic potential term).
5. **SYK check** (Sec. 2.5, Fig. 4): for N=2 SYK with q̂=3, N=11, 200 realizations, level statistics of the
   supercharge singular values in each multiplet sector match β=2, except the (−3/2, 3/2) sector, which is β=1
   (a particle–hole/CT effect).
6. N=4 JT: spectrum, two-boundary wormhole, beginnings of N=4 random matrix theory (Sec. 7).

## Relationship to our project
- Underlies CCSY's "supercharge chaos" conjecture and Johnson's model. It is the **energy-side** matrix model;
  it says nothing directly about BPS eigenvectors relative to UV operators such as Π_T.
- "Random BPS wavefunctions" in TW refers to randomness within the RMT model's own basis. The analogue for our
  decoder is BIY's Haar-random states *within* the BPS space, not Haar randomness of the BPS subspace inside the
  Fock space.
- The β=1 exception for particle–hole-symmetric sectors is a reminder to check the symmetry class of O_T at even N′,
  half filling. There, the audit found an exact λ ↔ 1−λ pairing. CV's ⟨r⟩ result is at N′=12, P=6, exactly such a
  sector; it still came out GUE, which should be understood.
- TR/loop-equation technology here (App. B) is the correct framework for any genuinely gravitational recursion,
  as opposed to TR applied to the Haar/Wachter curve in the draft's App. A.
