# Demni, Hamdi & Hmidi (2012) — On the spectral distribution of the free Jacobi process

**Citation.** N. Demni, T. Hamdi, T. Hmidi, *On the spectral distribution of the free Jacobi process*,
arXiv:1204.6227 [math.SP]. File: `papers/demni_hamdi_hmidi_2012.pdf`. (Read: Sec. 1 and statement of results.)

## Main content
- Free Jacobi process J_t = P Y_t Q Y_t* P with Y a free unitary Brownian motion, τ(P) = λθ, τ(Q) = θ.
  It interpolates from an initial (non-free) configuration at t = 0 to the free (Wachter) equilibrium as
  t → ∞ (Voiculescu's liberation process).
- Moment ODE (eq. 2): ∂_t m_n = −n m_n + θn m_{n−1} + λθn Σ_k m_{n−k−1}(m_k − m_{k+1}), m_n(t) = τ(J_t^n)/τ(P).
  In particular m₁(t) relaxes exponentially to θ.
- For λ = 1, θ = 1/2: explicit moments (eq. 3) and J_t ~ (1/4)(2 + Y_{2t} + Y_{2t}*); the support fills (0,1) at t = 2.

## Relationship to our project
Not yet used anywhere in the repo, but it gives the natural mathematical language for a *partially liberated*
pair of projections. Hypothesis (see `docs/research_questions.md`, Q2): the SYK decoder is a partial
liberation of the LES reference subspace B⁰ = B^p_N ⊕ (B^{p−1}_N ∧ e), whose decoder spectrum is exactly
Bernoulli on {0,1}, towards the Wachter law. The wall pile-up would then be the un-liberated remainder.
Caveat: free unitary Brownian motion is itself a Haar-type model. Whether SYK's rotation is "Brownian" is an open question.
