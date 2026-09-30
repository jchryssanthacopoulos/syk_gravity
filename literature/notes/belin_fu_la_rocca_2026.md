# Belin, Fu & La Rocca (2026) — Microcanonical operator spectrum, BPS chaos and free probability

**Citation.** A. Belin, Y. Fu, F. La Rocca, *Microcanonical Operator Spectrum, BPS Chaos and Free Probability
Theory*, arXiv:2609.22422v1 [hep-th], 18 Sep 2026 (Milano-Bicocca / GIST / PKU). File:
`papers/belin_fu_la_rocca_2026.pdf` (69 pp). Read in full on 2026-09-30, except the proofs in App. B–C and the
spin-1/4-level examples in App. D–E, which were skimmed.

## Main question
What is the eigenvalue distribution of a simple operator O once it is projected to a microcanonical window, or to a
degenerate (BPS) subspace, O_K = P_K O P_K? The paper assumes the window is Haar-randomly oriented relative to O's
eigenbasis. The motivation is BPS chaos (LMRS; Chen–Lin–Shenker), where the projected spectrum is "the only physical
quantity".

## Physical system
There is none specifically. The setting is an abstract N-dimensional Hilbert space, a fixed Hermitian O with given
spectrum f(O), and a rank-K projector P_K. Numerics use N = 1000.

## Important definitions
- O_K = P_K U O U† P_K with U Haar on U(N) (eqs. 1.2, 2.6). α = K/N (eq. 2.8).
- Ō = Tr O/N, σ_O² = Tr O²/N − Ō² (eqs. 2.9–2.10).
- Free compression: (pAp, φ/φ(p)) with p free from A (eqs. 4.66–4.67).

## Main assumptions
The Hamiltonian eigenbasis is Haar-random relative to O's eigenbasis (eq. 2.2): "chaos ⇒ random orientation". The
spectrum of O is fixed or i.i.d. and otherwise arbitrary, with mild scaling assumptions on its moments.

## Main analytical results
1. **Exact finite-N joint density.** This uses the HCIZ integral and a confluent Vandermonde limit. The density is a
   Vandermonde factor times det of cardinal B-splines M_{N−K+1}(O_j,…,O_{j+N−K}; λ_l) (eq. 2.33). The support
   condition O_j ≤ λ_j < O_{j+N−K} (eqs. 2.34–2.35) strengthens Cauchy interlacing.
2. **GUE central limit (their main result).** For N, K → ∞ with α ≪ 1, μ(λ) is GUE centred at Ō with variance set by
   σ_O² (eq. 2.65), up to O(α³). The one-point density is a semicircle with variance α σ_O² (eq. 2.66). Higher
   cumulants are suppressed as in eqs. 2.53–2.55.
3. **Sine kernel and level repulsion at finite α.** This is a heuristic argument (eqs. 2.67–2.77), valid when the
   projected density is smooth.
4. **Spin-½ operator (O = ±1, each with multiplicity N/2).**
   - For α ≤ ½, ρ(λ) = √(4α(1−α) − λ²)/(2πα(1−λ²)) (eq. 3.3, attributed to Iniguez–Srednicki and Collins). It is
     the arcsine law at α = ½ (eq. 3.7).
   - For α > ½ there is **eigenvalue recovery**: ±1 appear with multiplicity K − N/2 each (eq. 3.8).
5. **App. A, Theorems A.1–A.2.** An eigenvalue of O with multiplicity γ is an eigenvalue of O_K with multiplicity
   ≥ K − (N − γ). For a random (Haar or generic GL) rotation, equality holds almost surely: exceeding the bound is a
   measure-zero event.
6. **Free probability, §4.**
   - Thm. 4.1: A and UBU† are asymptotically free.
   - ρ_{O_K} = μ_{P} ⊠ f (eq. 4.33).
   - Thm. 4.2: moments agree with the semicircle to O(α²) and differ at O(α³) (eq. 4.36).
   - Second-order freeness: the covariance moments m_{p,q} agree with GUE to O(α) (eq. 4.65), checked for p, q ≤ 3.
   - Free compression μ_{pap} = μ_{αa}^{⊞1/α} (eq. 4.69) reproduces eq. 3.3 via its Cauchy transform (eq. 4.76).

## Numerical methods and results
Haar sampling with N = 1000 and 1000 samples: spin-½ and spin-1 operators, and a uniform random spectrum on [0, 1]
(Figs. 1–4, Table 1).
- **Degenerate O.** RMT statistics break down after eigenvalue recovery sets in.
- **Non-degenerate O.** The spacing distribution and SFF look GUE even for α close to 1.

## Discussion points relevant to BPS chaos (§5)
- They note (citing LMRS 2207.00408) that the simplest operators in SUSY SYK do **not** show the semicircle/GUE
  expected from Haar orientation.
- They say few-fermion operators are **not mutually free** after projection, citing *S. Khamnei & K. Papadodimas,
  unpublished* (their ref. [53]). They add that freeness improves for heavier operators, consistent with LMRS's super-JT.
- Open question they pose: "quantify how close to a true Haar-random unitary the physical rotation is", e.g. via
  k-designs (Roberts–Yoshida).

## Limitations
- Everything rests on the Haar assumption, so no specific Hamiltonian is studied. They explicitly call it a toy model.
- The free-probability proof covers only the one-point density and a few second-order moments.

## Relationship to our project
1. **It is our Haar null model, written from the BPS-chaos side.** With O = 2Π_T − 1 (spin-½, b = ½) and P_K = P_B
   (α = a):
   - Their eq. 3.3 is the Wachter law in our variable μ = 2λ − 1 at b = ½.
   - Their eigenvalue recovery (eq. 3.8) gives the Wachter wall atoms.
   - Their free-compression formula (eq. 4.69) is our free baseline (`nullmodels.free_compression_mu_moments`).
   - Their "BPS chaos ⇒ Haar orientation" hypothesis predicts r = a, the value our data exclude (r/a = 1.5–1.7 and
     growing; `docs/research_questions.md` Q1).
2. **SYK violates their generic-position theorem at the walls.** By Thms. A.1–A.2, a Haar-rotated BPS space has
   exactly max(0, d − D/2) eigenvalues pinned at each wall (b = ½), i.e. a wall-atom fraction (2a − 1)/a. The SYK
   decoder has more (`results/data/parity_family_2026-09-29/parity_family_table.md`, "atoms" column at s = 1):

   | N′ | 8 | 10 | 12 | 13 | 14 |
   |---|---|---|---|---|---|
   | SYK | 0.741 | 0.556 | 0.329 | 0.056 | 0.056 |
   | Haar (a.s.) | 0.704 | 0.444 | 0.099 | 0 | 0 |

   In their language, SYK's BPS space sits in a measure-zero, non-generic position relative to the slot. This is a
   clean, rigorous way to state our Q3 (atoms from structural rank bounds) [derived: Thm. A.2 plus our data].
3. **Our size spectrum answers their open question at second order** [derived]. Because the SYK coupling law is only
   U(N)-invariant, not U(D)-invariant, the second moment of the BPS-projector ensemble is fixed by the size
   spectrum: the U(N)-twirl of |P⟩⟩⟨⟨P| is Σ_k c_k Π_k with c_k = d w_k/dim V_k (multiplicity-free decomposition of
   End Λ^F). A U(D)-Haar projector, and more generally any U(D) 2-design, has flat c_k for k ≥ 1 in expectation.
   So:
   - w_k measures the **2-design defect** of the BPS-projector ensemble relative to Haar;
   - T₂ and the two-model fourth moment are second-moment quantities, fixed exactly by w_k;
   - the decoder's T₄ is a fourth-moment quantity, and its remainder R is where structure beyond the U(N)-twirl first
     enters.

   This is the quantitative answer to "how close to Haar is the physical rotation", for the most-studied BPS model.
4. **Regime.** Their GUE/semicircle result needs α ≪ 1. Our a = 0.4–0.8, so the Wachter law (finite α), not the
   semicircle, is the relevant Haar law at accessible N′. Asymptotically a → 0 and their regime applies to the
   *null model*, but SYK's r ≫ a shows the null model is not the right large-N description.
5. **Novelty check.** They cite unpublished work (Khamnei–Papadodimas) finding that few-fermion operators in SUSY SYK
   are not free after projection. That is qualitatively our "x-independent non-freeness" and "two BPS spaces are not
   free". Our exact finite-N characterization via size and transmission spectra goes further, but this is a
   **possible overlap to watch**.
6. Related references worth acquiring: Iniguez–Srednicki, *Microcanonical truncations of observables*
   (arXiv:2305.15702); Wang et al., *Emergence of unitary symmetry of microcanonically truncated operators*
   (arXiv:2310.20264).
