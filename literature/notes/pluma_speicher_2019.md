# Pluma & Speicher (2019) — A dynamical version of the SYK model and the q-Brownian motion

**Citation.** M. Pluma, R. Speicher, *A dynamical version of the SYK model and the q-Brownian motion*, Random
Matrices: Theory and Applications 11 (2022) 2250031, arXiv:1905.12999 (v3, 8 Sep 2020) [math.OA]. The journal
reference was taken from the reference lists of Boruch–Lin–Yan [11] and Lin 2022 [21]. File:
`papers/pluma_speicher_2019.pdf` (23 pp).

**Note written 2026-10-05.** Read in full.

## Main question
What is the large-n joint distribution of several independent SYK Hamiltonians? Does SYK with time-dependent
(Brownian) couplings converge to the q-Brownian motion? What are the fluctuations (higher cumulants of traces)?

## Physical system
Majorana SYK H_{n,q_n} = i^{⌊q_n/2⌋} binom(n,q_n)^{−1/2} Σ J_R Ψ_R (5). Here q_n is the interaction length (our p) and
q_n²/n → λ. Couplings are independent with mean 0, variance 1 and bounded moments; Gaussian for §4.

## Important definitions
- **q-Fock space.**
  - Bożejko–Speicher inner product; a(f)a*(g) − q a*(g)a(f) = ⟨f,g⟩.
  - q-Gaussians s_q(h) = a(h) + a*(h).
  - Vacuum state τ.
  - q-Wick formula: τ(s_q(f_{ε(1)})⋯s_q(f_{ε(k)})) = Σ_{π∈P₂(k)} q^{cr(π)} Π ⟨f_r, f_s⟩ (p. 7; eq. 9 for an
    orthonormal family).
- **q-Brownian motion:** S_q(t) = s_q(1_{[0,t]}) (10).

## Main assumptions
- q_n²/n → λ ∈ [0,∞].
- Parity of q_n is fixed: q = e^{−2λ} for even q_n, q = −e^{−2λ} for odd q_n (Thm. 3.1).
- Uniformly bounded coupling moments. Gaussian couplings for §4.

## Main analytical results
1. **Thm. 3.1 (multivariate).** p independent copies H₁,…,H_p converge in distribution to s_q(h₁),…,s_q(h_p) with
   h_i orthonormal: lim E tr(H_{ε(1)}⋯H_{ε(k)}) = Σ_{π ≤ ker ε} q^{cr(π)} (11).
   - Cross-copy crossings carry the same q as same-copy crossings.
2. **Cor. 3.2 (dynamical).** Couplings J_R(t) that are independent classical Brownian motions give H(t) → q-Brownian
   motion (13).
   - The proof is one line: expand in orthogonal increments and use **multilinearity**.
3. **Lemma 3.3.** tr(Ψ_α) = (−1)^{q_n cr(π) + Σ|α(V)∩α(W)|} for ker α ≥ π (14). This is the exact sign of a word;
   each crossing contributes (−1)^{q_n + overlap}.
4. **Lemma 3.4(ii).** A restricted sum with one block value fixed is independent of that value (permutation
   symmetry of [n]).
5. **Thm. 4.1 (fluctuations, higher cumulants).**
   - |I_n|^{m−1} c_m(tr H_{ε₁},…,tr H_{ε_m}) → Σ q^{cr(π′)} over pairs (π, π′) with π ∨ θ = 1, π′ ≤ θ,
     |π ∨ π′| = k/2 − m + 1 (22).
   - Sample-to-sample fluctuations of traces are therefore O(|I_n|^{−1/2}), with |I_n| = binom(n,q_n).
   - Remark 4.2 notes a structural resemblance to the partitioned permutations of higher-order freeness.
6. **§6.**
   - There is *no* good analytic description of the multivariate q-Gaussian law.
   - The operator-valued Cauchy transform has no operator-valued continued fraction (Anshelevich–Williams). This is
     cited as the obstacle to closed-form multi-flavour (e.g. BINT 4-pt) computations.

## Important equations for us
(9), (11), (13), (14), (22).

## Numerical methods
None.

## Limitations
- Majorana, non-SUSY.
- Moments of Hamiltonians only (polynomials). No spectral projectors, no ground states.
- Leading order in n.

## Relationship to our project
1. **Correlated couplings reduce to independent ones by multilinearity (Cor. 3.2 argument).**
   - Any jointly Gaussian family of coupling vectors with the same correlation C_ij for every index set gives, in
     the limit, s_q(h_i) with ⟨h_i, h_j⟩ = C_ij.
   - For OU-correlated copies ω_x = xω + √(1−x²)ω′ this means (H(ω), H(ω_x)) → (s_q(h), s_q(h_x)) with
     ⟨h, h_x⟩ = x. Each chord between the copies weighs x, and every crossing weighs q irrespective of species.
   - **Correction to our BBNR note.** That note called OU-correlated copies "a correlated mixed q-Gaussian, beyond
     Bożejko–Speicher and Liu–Shen". This is wrong for the Majorana case: it is the ordinary (mixed) q-Gaussian
     functor evaluated on non-orthogonal vectors. The note has been corrected.
   - For N = 2 (complex couplings, oriented chords) the same multilinearity argument should apply to BBNR's
     X/O-chord rules. We have not seen this written down. [plausible, not established]
2. **The OU process is the q-OU process.**
   - Brownian couplings give q-Brownian motion. Stationary OU couplings (covariance e^{−|t−s|}) give the q-Ornstein–
     Uhlenbeck process of Bożejko–Kümmerer–Speicher, which is Markov (BKS Cor. 3.10).
   - Our T₂^OU(x) is its two-time (time-ordered) correlator. T₄^OU(x) = E τ(P(ω)P(ω_x)P(ω)P(ω_x)) is **out of OU-time
     order** (s, t, s, t), so the Markov kernel alone does not determine it. This is the precise sense in which R
     belongs to a crossed/OTOC channel. See BKS note and `docs/derivations.md` D5.
3. **Fluctuation scale.**
   - Sample-to-sample fluctuations of normalized traces are suppressed by binom(N,p)^{−1/2} (Thm. 4.1). At N′ = 14,
     p = 3, binom = 364.
   - This sets the expected size of realization-to-realization variation in per-realization quantities such as the
     twirl term S. It is relevant to judging whether ⟨S⟩ vs per-realization S matters.
4. **§6 caution.** No closed-form analytic machinery exists even for two independent q-Gaussians. For T₄ we should
   expect to compute in the chord Fock space (transfer matrices) rather than through a closed-form R-transform.
