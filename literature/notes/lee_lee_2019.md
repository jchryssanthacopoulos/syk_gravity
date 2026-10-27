# Lee & Lee (2019) — Notes on superconformal representations in two dimensions

**Citation.** S. Lee, S. Lee, *Notes on superconformal representations in two dimensions*, Nucl. Phys. B 956 (2020) 115033,
arXiv:1911.10391v2 [hep-th] (9 Oct 2020). The journal reference is as given in G. Lin 2026 ref. [61], not on the PDF.
File: `papers/lee_lee_2019.pdf`.

What was read (2026-10-27):
- abstract and contents;
- §2.2 (multiplets, unitarity, notation);
- §2.3 (top components);
- §3.1–3.2 (N = 1, 2 tables);
- §4.1 (recombination rules, the generic and N ≤ 4 cases);
- §6.2.1 (N = 2 characters).

The higher-N sections were not read.

**Role in this project.** Representation-theory reference for the global superconformal algebras, chiral half. For
N = 2 the chiral global algebra (L₀, L_{±1}, G^±_{±1/2}, J₀) is su(1,1|1). This is the gauge superalgebra of N=2 super-JT
(Penington–Witten), so Lee–Lee's N = 2 tables classify the matter multiplets that appear in PW and in G. Lin 2026.

**Correction to our earlier synthesis.** `synthesis.md` §6 said this paper supplies "superconformal representation and
exchange data needed for a closed-form super-JT OTOC (G. Lin 2026 §6)". That overstates it.
- G. Lin 2026 §6 cites [61] for:
  - the *notation* A[2J]_h / L[2J]_h;
  - the *channel sum*, i.e. which supermultiplets appear in the OPE of two external multiplets.
- G. Lin's own example is small N = 4: A[1]_{1/2} ⊗ A[1]_{1/2} = A[2]₁ ⊕ ⊕_{m≥0} L[0]_{1+m}.
- Lee–Lee contain no tensor-product decompositions, OPE coefficients, 6j symbols or exchange matrices. A search of the text
  for "tensor product", "OPE", "6j", "crossing" and "fusion" finds only R-symmetry tensor products.
- Channel sums can be *derived* from their characters, by multiplying characters and decomposing. The exchange matrix Λ_Ξ
  that G. Lin needs for the OTOC (his eq. 6.8) is not here, and G. Lin says it is not yet known.

## Main question
Classify the unitary multiplets of the global subalgebras of the 2D superconformal algebras, for all N. Then:
- decompose them into conformal multiplets (Racah–Speiser);
- give recombination rules and characters;
- apply them to stress tensors, conserved currents, deformations and supersymmetry enhancement.

## Physical system
Representation theory of osp(N|2)-type global algebras: one chiral half of a 2D SCFT, with the right movers treated
identically.

## Important definitions (§2.2)
- Superconformal primary [R]_h: R-symmetry label [R] and L₀-weight h.
- Unitarity bound from the first-level states (2.10)–(2.11).
- Long multiplet L[R]_h (h > h₁(R)); short multiplet A[R]_h (h = h₁(R)); vacuum V[0]₀.
- Components are conformal primaries obtained by acting with G_{−1/2}. Level l has weight h₀ + l/2.
  - Long multiplets end at level N; short ones end earlier.
- Generic and sporadic top components are the supersymmetric (all-G-annihilated) components (§2.3).

## Main results used here
1. **N = 2 (§3.2, Table 2, eq. 3.2).** The U(1) R-charge j₀ has G^± of charge ±1. The unitarity bound is h₀ ≥ |j₀|/2.
   - Long: L[j₀]_{h₀}, h₀ > |j₀|/2, with generic top [j₀]_{h₀+1}.
   - Short, j₀ > 0: A[j₀]_{j₀/2}. Null component [j₀+1]_{j₀/2+½}; generic top [j₀−1]_{j₀/2+½}. Mirror for j₀ < 0.
   - Vacuum: V[0]₀, with null components [±1]_{½}.
   - The global bound is weaker than the super-Virasoro one (text after 3.2).
2. **Recombination (§4.1).** In 2D every short multiplet is a limit of a long one (end of §2.2), so no multiplet is
   absolutely protected. The N = 2 rule (4.4):
   - L[j₀>0]_h → A[j₀]_{j₀/2} ⊕ A[j₀+1]_{j₀/2+½} as h → j₀/2;
   - L[0]_h → A[1]_{½} ⊕ A[−1]_{½} as h → 0.
   - Generic rule (4.2); small N = 4 (4.5).
   - Consistency check: the state count matches on both sides.
3. **Characters (§4.2; N = 2 in §6.2.1).**
   - Massive: q^{h−c/24} y^j (1+q^{½}y)(1+q^{½}y^{−1})/(1−q) (6.10).
   - Massless: the same times 1/(1+q^{½}y^{sign j}) (6.11).
   - The recombination rule in character form is (6.12).

## Dictionary to Penington–Witten (derived here; not stated in either paper)
PW label an SU(1,1|1) multiplet by (λ, u): λ is the lowest j₁ eigenvalue and u the R-charge of that state.
- In PW's BPS representation (u > 0) the R-charges are j = u and u − ½ (4.116).
- In Lee–Lee's A[j₀]_{j₀/2} they are j₀ and j₀ − 1.

So **λ = h₀ and u = j₀/2**, with PW's long condition |u| < λ ⇔ Lee–Lee's h₀ > |j₀|/2.

The PW zero-energy condition |2u − J_R| < ½ then reads |j₀ − J_R| < ½ in Lee–Lee units. This matching uses only the
charge assignments of the components and could fail if a sign convention for the R-current differs; the content
(long ↔ four SL(2) series, BPS ↔ two) agrees.

## Numerical methods
None.

## Limitations
- Global multiplets only. The super-Virasoro decomposition is given separately (§6), but the 1D/JT problem needs only the
  global part.
- The paper conjectures that first-level unitarity suffices (end of §2.2).
- No OPE, fusion or exchange data (see the correction above).

## Relationship to our project
- **Roadmap C3/A3 (an N=2 OTOC via G. Lin's OPE route).** Lee–Lee provide the multiplet list and characters, from which
  the N=2 channel sums can be computed. The remaining hard input, the super-6j/exchange matrix, is not in this paper. If
  we pursue the OTOC route, that input must be derived (or found elsewhere).
- **Recombination as an algebraic picture of fragility.**
  - The rule L → A ⊕ A′ says that two short multiplets whose charges differ by one can pair up into a long multiplet
    under a deformation.
  - In N=2 SYK the exceptional BPS states (CV 2026; D10) come in index-cancelling pairs of opposite fermion parity.
    CV's q = 3 example (§2.2) has them at degrees p = 1 and p = 4.
  - [speculative analogy] The fragile/exceptional pairs are the QM analogue of A ⊕ A′ pairs that may recombine. Monotone
    states would then correspond to shorts with no partner. This is only an analogy: SYK Q-multiplets are two-state
    {Q, Q†} doublets, not su(1,1|1) multiplets, and no superconformal symmetry acts on the microscopic Hilbert space.
    Worth testing only if a chord or bulk construction assigns su(1,1|1) labels to uplift states (roadmap C).
