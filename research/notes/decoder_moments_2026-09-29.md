# Higher moments of the decoder: parity form, correlated BPS spaces, and an exact size sum rule

*2026-09-29. Foothold note for open problem 2 of the checkpoint (`research/pdfs/checkpoint_2026-09-28.pdf` §8):
"higher moments and the density". Status labels: **[derived]**, **[verified]** (exact numerics confirm a derived
statement), **[num]** (numerical observation), **[conj]**, **[falsified]**.*

Code: `src/parity_decoder.py`, `src/size_decomposition.py`, `src/nullmodels.py`; drivers
`scripts/run_parity_family.py`, `scripts/run_two_model_overlap.py`, `scripts/run_size_decomposition.py`; analysis
`scripts/analyze_parity_family.py`, `scripts/analyze_size_sumrule.py`; tests `tests/test_size_decomposition.py`.
Data: `results/data/parity_family_2026-09-29/`. Saved tables: `parity_family_table.md`, `size_sumrule_tables.md`.
All runs under `scripts/memwatch.py` (peaks ≤ 0.9 GB).

## 0. Summary

1. **[derived] Parity form.** The decoder is $(1+PUP)/2$ with $U=(-1)^{n_i}$. Since $U^2=1$, $(PUP)^2 = PP'P$ with
   $P' = UPU$. So the decoder law is the principal-angle law between the BPS spaces of two SYK supercharges,
   $Q$ and $Q'=UQU$. $Q'$ is $Q$ with the sign of every coupling containing site $i$ flipped, so in chord language
   the two supercharges are *correlated* with coefficient $x$ per chord.
2. **[derived] Z₂ chord rule.** In a trace with $m$ insertions of $U$, a Q-chord gets weight $x$ iff it separates
   an odd number of $U$'s, and weight 1 otherwise. For $m=2$ this is BLY's horizontal matter chord. For $m\ge4$ it
   is *not* the Gaussian-matter rule.
3. **[num] Tunable-$x$ family.** Parity projectors on $s$ sites give $x_s$ from $\approx0.57$ down to exactly 0 at
   $s=N'/2$, with the same BPS projector.
   - As $x\to0$: $r\to a$ (the Haar value), and the second moment becomes free.
   - **But $T_4,T_6$ stay above their free values at $x=0$** ($+5.6\%$, $+10\%$ at $N'=14$, growing with $N'$).
4. **[num] Two independent SYK models** have BPS spaces that are *not* mutually free. $E[\nu]=a$ exactly (proved,
   §5), but $E[\nu^2]$ exceeds the free value by $0.6, 2.0, 4.3\%$ at $N'=10, 12, 14$ (16, 4, 16 pairs; exact
   prediction in §10). This is the $x$-independent part
   of the decoder's non-freeness.
5. **[derived + verified] Exact finite-$N$ sum rule** (§6). The coupling ensemble is $U(N)$-invariant. Decompose
   the BPS projector by operator size, $P=\sum_kP^{(k)}$ ($U(N)$-irreps $V_k\subset\mathrm{End}(\Lambda^P)$). Then
   $$E[T_2]=\sum_{k}w_k\,\hat\chi_k(U),\qquad w_0=a\ \text{exactly},$$
   with $\hat\chi_k$ the normalized $U(N)$ character, known in closed form. This is an **exact finite-$N$ version of
   the wormhole identity** (checkpoint eq. 15):
   - the Haar/free value is the size-0 (identity) component of $P$;
   - the deviation from freeness is $P$'s weight on small nonzero sizes.

   It reproduces the campaign's $r$ to 0.1–0.8 % for $N'=8\ldots14$, including odd $N'$. The size weights $w_{2n}$
   track BLY's chord-length distribution $P_n$, and $\hat\chi_{2n}$ tracks $x^{2n}$.
6. **[derived + verified] Exact two-copy structure (§10, added 2026-09-29).** The *transmission spectrum*
   $\phi_k$ (the eigenvalue of $X\mapsto E[PXP]$ on $V_k$) is a universal linear function of the size spectrum
   $w_k$. The fourth moment of two independent models is then exact:
   $$E\,\tau(P_1P_2P_1P_2)=\sum_k\phi_kw_k=(2a^2-a^3)+\sum_{k\ge1}(\phi_k-a^2)\,w_k .$$
   It matches independent-pair measurements to 0.01–0.2 % ($N'=10$–14). Haar projectors have flat transmission
   $\phi_{k\ge1}=a^2$. SYK has excess transmission $\delta\phi_k=\phi_k-a^2>0$ at small sizes, and the decoder
   variance is exactly $r=\phi_1/a$ at half filling. **Both kinds of non-freeness are controlled by $\delta\phi$.**
7. **[falsified, strong form] "Wachter = zero-length truncation of the chord computation for all moments"**
   (checkpoint §7.4 item 3). Removing the matter dressing ($x\to0$) does *not* leave the free law beyond $m_2$.
   **[conj, refined]** Free compression requires both $x=0$ *and* Haar-like multi-copy size structure of the BPS
   projector. SYK has neither.

## 1. Exact reformulation [derived]

$\Pi=1-n_i=\tfrac12(1+U)$ with $U=(-1)^{n_i}$, a unitary involution. So $P\Pi P=\tfrac12(P+PUP)$, and on $B$ the
eigenvalues satisfy $\lambda=(1+\mu)/2$, $\mu\in\mathrm{spec}(PUP|_B)\subset[-1,1]$. Define
$T_m=\tau\big((PUP)^m\big)=E[\mu^m]$. Then the decoder moments are
$$m_k=E[\lambda^k]=2^{-k}\sum_{m=0}^k\binom km T_m .$$
- At half filling (even $N'$), $T_{\rm odd}=0$ by particle–hole symmetry (§6.4).
- $r=\mathrm{Var}(\lambda)/[b(1-b)]=(T_2-c^2)/(1-c^2)$ with $c=2b-1$ ($c=0$ at even $N'$, $c=1/N'$ at odd $N'$).

**Correlated BPS spaces.** $U^2=1$ gives $(PUP)^2=P(UPU)P=PP'P$, where $P'=UPU$ is the BPS projector of
$Q'=UQU$. Writing $Q=\sum_IC_I\psi_I$,
$$Q'=\sum_I(-1)^{[i\in I]}C_I\psi_I .$$
Hence $\mu^2$ is the squared cosine of the principal angles between $B=\ker Q\cap\ker Q^\dagger$ and $B'=UB$.
Averaged over the random index sets, the per-term coupling correlation between $Q$ and $Q'$ is
$x=E[(-1)^{[i\in I]}]$.

## 2. The Z₂ chord rule [derived]

For a supercharge term $\psi_I$ ($|I|=p$), $U\psi_IU=(-1)^{[i\in I]}\psi_I$. Consider a trace of a word in
$Q,Q^\dagger$ with $m$ insertions of $U$. Moving all $U$'s together leaves $\mathrm{tr}(U^m)$ times a sign
$(-1)^{[i\in I]}$ for each supercharge factor, once for every $U$ it is moved past.
- A Q-chord (a Wick pair $C_IC_I^*$) therefore picks up $(-1)^{[i\in I]\cdot\#}$, where $\#$ is the number of $U$'s
  its two endpoints are separated by.
- Averaging over $I$ gives weight **$x$ if $\#$ is odd, 1 if even**:
  - exact finite-$N$: $x=1-2p/N$;
  - double-scaled: $x=e^{-2p/N}=q^{1/p}$, the BLY matter weight $q^\Delta$ with $\Delta=1/p$.
- Equivalently, label the $m$ zero-temperature arcs alternately $P$ and $P'$. Chords joining a $P$-arc to a
  $P'$-arc get $x$; chords joining arcs of the same kind get 1. This is the chord theory of two supercharges with
  mixed-contraction weight $x$.
- $m=2$: chords crossing between the two arcs get $x$ each, so $T_2=\sum_nP_nx^{2n}$ (checkpoint eq. 15). ✓
- $m=4$: a chord from arc 1 to arc 3 (same kind) gets 1, not $x^2$. So $T_{2k\ge4}$ differs from the random-matter
  (q-Gaussian) answer, and from a naive sum over pairings.

## 3. The parity family: tunable $x$ with the same $P$ [num]

Replace $U$ by $U_S=(-1)^{N_S}$ with $|S|=s$ and $S\ni N'-1$. Then $U_S\psi_IU_S=(-1)^{|I\cap S|}\psi_I$ and
$$x_s=\sum_k(-1)^k\binom sk\binom{N-s}{p-k}\Big/\binom Np\quad\text{(exact)},\qquad x_s\approx e^{-2ps/N}\quad(s\ll N).$$
For odd $p$, even $N$ and $s=N/2$, $x=0$ **exactly** (complement symmetry: $U_{S^c}=(-1)^{N_F}U_S$ with $N_F$ fixed).
The BPS projector is unchanged, so this isolates the role of $x$.

$N'=14$, $q=3$, 4 SYK realizations; Haar controls give $1.000\pm0.001$ in every column
(`parity_family_table.md`):

| $s$ | $x_s$ exact | $x_s$ double-scaled | $r/a$ | $T_4/T_4^{\rm free}$ | $T_6/T_6^{\rm free}$ | KS to Wachter |
|---|---|---|---|---|---|---|
| 1 (decoder) | +0.571 | 0.651 | 1.518 | 1.887 ± 0.018 | 2.21 | 0.145 |
| 2 | +0.275 | 0.424 | 1.208 | 1.359 ± 0.024 | 1.49 | 0.080 |
| 3 | +0.088 | 0.276 | 1.062 | 1.144 ± 0.024 | 1.22 | 0.035 |
| 4 | −0.011 | 0.180 | 1.007 | 1.056 ± 0.014 | 1.11 | 0.026 |
| 5 | −0.044 | 0.117 | 1.004 | 1.049 ± 0.006 | 1.10 | 0.021 |
| 7 | 0 | 0.050 | 1.010 | 1.056 ± 0.002 | 1.10 | 0.021 |

At $x\approx0$ (last row), across sizes:

| $N'$ | $r/a$ | $T_4/T_4^{\rm free}$ | $T_6/T_6^{\rm free}$ | KS (SYK / Haar) |
|---|---|---|---|---|
| 10 | 1.010 | 1.014 ± 0.007 | 1.016 | 0.014 / 0.012 |
| 12 | 1.004 | 1.025 ± 0.004 | 1.043 | 0.018 / 0.005 |
| 13 | 0.997 | 1.032 ± 0.010 | 1.069 | 0.019 / 0.004 |
| 14 | 1.010 | 1.056 ± 0.002 | 1.102 | 0.021 / 0.002 |

($T_4^{\rm free}=2a^2-a^3$ and $T_6^{\rm free}=5a^3-6a^4+2a^5$ at $b=\tfrac12$: free compression of a symmetric
$\pm1$ variable, whose free cumulants are $\kappa_2=1$, $\kappa_4=-1$, $\kappa_6=2$.)

Findings:
- **The second moment** becomes free as $x\to0$ ($r/a\to1$), as the zero-length identity requires. The small
  residual $+0.4\ldots1\%$ is explained exactly in §6.3.
- **Higher moments do not.** An $x$-independent excess remains, and it grows with $N'$.
- **Decomposition of the decoder's $T_4$ excess at $N'=14$:** $T_4^{\rm dec}-T_4^{\rm free}=0.252$, of which
  $0.236$ (94 %) disappears as $x\to0$ and $0.016$ (6 %) does not.
- **Extracting $P_n$ from $T_2(x_s)$ is ill-posed.** With exact $x_s$ it gives $P_1>1$. With double-scaled $x_s$ the
  intermediate-$s$ prediction is 5–10 % high. The effective per-chord weight at finite $N$ lies between the two;
  §6 removes this ambiguity.

## 4. Two independent SYK models [num]

Take two independent coupling draws (same $N'$, same $P$), and let $\nu$ be the squared cosines of the principal
angles between $B_1$ and $B_2$, so $E[\nu^m]=\tau((P_1P_2P_1)^m)$. Free compression gives
$\nu\sim\mathrm{Wachter}(a,a)$. With 4 seed pairs per $N'$ (`size_sumrule_tables.md` §4):

| $N'$ | $E[\nu]/a$ | $E[\nu^2]/{\rm free}$ | $E[\nu^3]/{\rm free}$ | KS to Wachter(a,a) (SYK / Haar) |
|---|---|---|---|---|
| 10 | 1.009 ± 0.006 | 1.014 ± 0.005 → **1.006 (16 pairs, CV)** | 1.017 ± 0.004 | 0.026 / 0.010 |
| 12 | 0.999 ± 0.003 | 1.021 ± 0.004 | 1.040 ± 0.004 | 0.032 / 0.007 |
| 14 | 0.993 ± 0.002 | 1.034 ± 0.003 → **1.043 (16 pairs, CV)** | 1.077 ± 0.004 | 0.038 / 0.007 |

- **Independent SYK BPS spaces are not free relative to each other,** and the deviation grows with $N'$.
- **They are close to the $x=0$ parity numbers** ($T_4$ ratio: parity 1.014 / 1.025 / 1.056 vs two-model
  1.006 / 1.020 / 1.043 at $N'=10/12/14$). The parity decoder keeps a small extra residual correlation (the
  nonzero $\hat\chi_{k\ge2}(U_{N/2})$ of §6.3).
- The bold entries (16 pairs, control variate on the exact $E[\nu]=a$) supersede the 4-pair values; see §10.5.
- **[interp] Mechanism.** In chord language the two models' chords still cross each other with the $q$-weights set
  by shared site indices, so the models are not free unless $q\to0$. The growth with $N'$ (smaller $\lambda$, larger
  $q=e^{-\lambda}$) matches this. Representation-theoretically (§6), both projectors concentrate on low operator
  sizes, which a Haar projector does not.

## 5. $U(N)$ invariance and $E[\nu]=a$ [derived]

The couplings $C_{ijk}$ are i.i.d. standard complex Gaussians in the orthonormal basis $e_i\wedge e_j\wedge e_k$
of $\Lambda^3\mathbb C^N$. The unitary action of $U(N)$ on single-particle modes induces a unitary action on
$\Lambda^3$, so the coupling law is $U(N)$-invariant. $Q=\omega\wedge$ is covariant, so $P_B(g\omega)=gP_B(\omega)g^\dagger$,
and the BPS-projector ensemble is invariant under conjugation by $U(N)$ (acting on $\Lambda^P$). Since $\Lambda^P$
is irreducible, Schur's lemma gives
$$E[P]=a\,\mathbb 1\qquad\text{exactly}.$$
Hence $E_2[\tau(P_1P_2)]=\tau(P_1\,a\mathbb 1)=a$ for every $P_1$. The data ($0.993$–$1.009$) are consistent at
$\lesssim3\sigma$ with 4 pairs.

## 6. Exact finite-$N$ sum rule: operator size of the BPS projector [derived + verified]

### 6.1 Setup
$\mathrm{End}(\Lambda^P)=\Lambda^P\otimes(\Lambda^P)^*$ decomposes under $U(N)$ multiplicity-free as
$\bigoplus_{k=0}^{\min(P,N-P)}V_k$. Here $V_k$ has highest weight $(1^k,0,\ldots,0,(-1)^k)$ ("traceless $k$-body
operators") and $\dim V_k=\binom Nk^2-\binom N{k-1}^2$. The adjoint Casimir
$C(X)=\sum_{ij}[E_{ij},[E_{ji},X]]$ ($E_{ij}=\bar\psi_i\psi_j$) acts on $V_k$ as $2k(N+1-k)$. The normalized character
of $V_k$ at $g\in U(N)$ is $\hat\chi_k(g)=(|\chi_{\Lambda^k}(g)|^2-|\chi_{\Lambda^{k-1}}(g)|^2)/\dim V_k$. For the parity
$U_S$ with $|S|=s$, $\chi_{\Lambda^k}(U_S)$ is the Krawtchouk polynomial $K_k(s)=\sum_j(-1)^j\binom sj\binom{N-s}{k-j}$.

**Verified** (test `test_casimir_spectrum_and_characters`): full diagonalization of the superoperator at $N=6,7$
gives eigenvalues $2k(N+1-k)$ with multiplicities $\dim V_k$. $\mathrm{Tr}_{V_k}\mathrm{Ad}_{U_S}/\dim V_k$ equals the
Krawtchouk formula exactly for $s=1,2$.

### 6.2 The sum rule
Decompose $P=\sum_kP^{(k)}$ with $P^{(k)}\in V_k$ (Frobenius-orthogonal), and set $w_k=\|P^{(k)}\|^2/d$. Then:
1. $\sum_kw_k=\mathrm{Tr}P^2/d=1$.
2. **$w_0=a$ exactly, for every realization.** $P^{(0)}$ is the $U(N)$-invariant part of $P$, which is
   $(\mathrm{Tr}P/D)\mathbb 1=a\mathbb 1$, with norm$^2$ $a^2D=ad$.
3. $U_S\in U(N)$, so $\mathrm{Ad}_{U_S}$ preserves every $V_k$, and
   $T_2=\tau(PU_SPU_S)=\frac1d\sum_k\langle P^{(k)},\mathrm{Ad}_{U_S}P^{(k)}\rangle$ with no cross terms.
4. By $U(N)$ invariance, $E\big[|P^{(k)}\rangle\rangle\langle\langle P^{(k)}|\big]$ commutes with $U(N)$ on the
   irreducible $V_k$, so by Schur it is proportional to $\mathbb 1_{V_k}$. Hence
$$\boxed{\;E[T_2(S)]=\sum_{k\ge0}\bar w_k\,\hat\chi_k(U_S)=a+\sum_{k\ge1}\bar w_k\,\hat\chi_k(U_S)\;}$$

This is exact at finite $N$, for any $S$ and any charge sector, with no double scaling and no free parameter.

**Haar check [derived].** For Haar $P$, $\bar w_k=(1-a)\dim V_k/(D^2-1)$ for $k\ge1$ (observed: e.g. $w_5=0.2127$ vs
$0.212$ at $N'=12$). Using $\sum_k\dim V_k\hat\chi_k=|\mathrm{Tr}\,U|^2$, at half filling
$E[T_2]=(aD^2-1)/(D^2-1)$, which is exactly the Weingarten $m_2$ of the audit converted to $\mu$
(test `test_haar_sum_rule_reproduces_weingarten`).

### 6.3 Results [verified]
The weights $w_k$ are computed with Lagrange filters in the Casimir; the residual
$\max|\sum_kP^{(k)}-P|\le2\times10^{-15}$.

| $N'$ | $w_0$ (= $a$) | $w_2$ | $w_4$ | $w_6$ | odd $k$ | $\sum_kw_k\hat\chi_k(U)$ | campaign $r$ |
|---|---|---|---|---|---|---|---|
| 8 | 0.7714 | 0.2128 | 0.0158 | – | 0 | 0.8170 | 0.8186 ± 0.0008 |
| 10 | 0.6429 | 0.2997 | 0.0574 | – | 0 | 0.7454 | 0.7463 ± 0.0007 |
| 12 | 0.5260 | 0.3486 | 0.1119 | 0.0135 | 0 | 0.6855 | 0.6874 ± 0.0009 |
| 13 | 0.4248 | 0.3398 | 0.1435 | 0.0318 | 0.015, 0.022, 0.024 | 0.6176 | 0.6171 ± 0.0018 |
| 14 (4 seeds) | 0.4248 | 0.3706 | 0.1629 | 0.0416 | 0 | 0.6371 | 0.6420 ± 0.0037 |

- **Single-site decoder.** The parameter-free sum rule reproduces $r$ to 0.1–0.8 %. Per realization with one
  fixed site, $T_2$ scatters by $\pm0.02$ around the prediction. Averaging over the 6 campaign modes of seed 0 at
  $N'=14$ gives $0.6375\pm0.008$ against $0.6370$.
- **Parity family.** Averaged over realizations it follows the same rule for every $s$. At $s=N'/2$ the per-chord
  weight is $x=0$, but $\hat\chi_{k\ge2}(U_{N/2})$ is small and nonzero ($K_k(N/2)=0$ for odd $k$,
  $\pm\binom{N/2}{k/2}$ for even $k$). That predicts $r/a=1.006$ at $N'=14$, against 1.010 measured, and explains the
  residual in §3.
- **Haar projectors** put almost all weight at large $k$, where $\hat\chi_k\approx0$, so $E[T_2]=a$ up to
  $O(D^{-2})$.

### 6.4 Only even sizes at half filling [derived (sketch) + verified]
- At even $N'$ ($P=N'/2$) every odd-$k$ weight vanishes (to $10^{-4}$ printed; exactly in the derivation). At odd
  $N'$ it does not (table).
- Derivation (sketch; signs not tracked in detail): particle–hole conjugation $V$ ($\psi\leftrightarrow\bar\psi$) maps $\Lambda^{N/2}$ to itself and
  $Q=\omega\wedge$ to the supercharge of the conjugate model, so $VB(\omega)=B(\bar\omega)=\overline{B(\omega)}$, i.e.
  $VPV^\dagger=P^T$. On a traceless $k$-body operator, $X\mapsto (VXV^\dagger)^T$ acts as $(-1)^k$. Therefore
  $P^{(k)}=0$ for odd $k$.
- The same map sends $n_i\to1-n_i$, which is the $\lambda\leftrightarrow1-\lambda$ pairing noted in the audit.

### 6.5 Relation to the chord theory [num + interp]
- **The size distribution tracks BLY's wormhole-length distribution.** Chord number $n$ corresponds to size
  $k=2n$ (`size_sumrule_tables.md` §3). At $N'=14$: $w_{2n}=(0.425,0.371,0.163,0.042)$ vs
  $P_n=(0.459,0.351,0.134,0.041)$.
- **The characters play the role of the matter weights:** $\hat\chi_{2n}(U)=(1,0.505,0.162,-0.029)$ vs
  $x^{2n}=(1,0.424,0.180,0.076)$.
- **Termwise errors of 10–20 % partly cancel in the sum,** which is why BLY's double-scaled formula is good to 2 %.
- **[interp] This is the exact finite-$N$ content of checkpoint eq. (15).**
  - The zero-length wormhole is the identity component of the BPS projector. It is exactly $a\mathbb 1$, which
    resolves the $a_\lambda$-vs-$a$ caveat.
  - A wormhole of length $n$ is the size-$2n$ component.
  - The matter propagator $q^{2\Delta n}$ is the normalized character of the probe's group element on that sector.
  - **Freeness (Haar) means the BPS projector's non-identity weight is spread uniformly over all operator sizes
    ($\propto\dim V_k$).** SYK's weight is concentrated at small sizes, which a local probe cannot randomize.
- The identification $n\leftrightarrow k/2$ is an empirical match so far. It is not derived.

## 7. What this says about the project's framing

- **Second moment: settled at finite $N$.** The deviation from freeness of the second free cumulant *is* the
  small-size content of the BPS projector, by an exact theorem (§6.2). The chord/wormhole reading is its
  double-scaled image.
- **Higher moments.** The decoder's non-freeness has two parts:
  - (i) $x$-dependent dressing: 94 % of the $T_4$ excess at $N'=14$;
  - (ii) an $x$-independent part: the BPS spaces of *different* SYK supercharges are not mutually free (§4), and
    this part grows with $N'$.
- **So "free compression = zero-length truncation" holds for $m_2$ only [falsified in strong form].**
- **[conj, refined]** The full free law requires Haar-like size structure of the *multi-copy* projector moments
  ($E[P^{\otimes2}]$, …), not just of $P$. Part (ii) should be governed by the two-copy analogue of $w_k$.

## 8. Next steps (in order)

1. **Exact $T_4$ from two-copy size data — DONE for two independent models (§10).** Original plan, kept for the record: For two independent models,
   $E\,\mathrm{Tr}(P_1P_2P_1P_2)=\mathrm{Tr}\big[E(P^{\otimes2})\,E(P^{\otimes2})\,\mathrm{SWAP}\big]$.
   - $E(P^{\otimes 2})$ lies in the $U(N)$ commutant on $\Lambda^P\otimes\Lambda^P$, which is multiplicity-free,
     indexed by $k=0\ldots P$.
   - So $E\tau(P_1P_2P_1P_2)=\frac1d\sum_\kappa\varepsilon_\kappa\,(E\,\mathrm{Tr}[P^{\otimes2}\Pi_\kappa])^2/\dim\kappa$
     ($\varepsilon_\kappa=\pm1$, the SWAP eigenvalue). The ingredients are computable from single-model data.
   - Test against §4. Then generalize to $\tau((PUP)^4)$, which needs $E[P^{\otimes2}]$ against $U^{\otimes2}$-twisted
     operators on $\mathrm{End}(\Lambda^P)^{\otimes2}$.
2. **Chord side.** The correlated two-supercharge chord Hilbert space (BLY §5 machinery with a Z₂ flag per chord,
   §2 above) gives $T_4(x)$ at $x=0$ and general $x$. Compare with 1 and with the data.
3. **Derive $n\leftrightarrow k/2$.** Relate chord number to operator size for $\mathcal N=2$ (BLY §5.2 relate
   $n_O,n_X$ to operator sizes), including why each $XO$ pair carries size 2.
4. **Density.** Once $T_2,T_4,T_6$ have exact/sum-rule descriptions, reconstruct the density. Compare with Wachter,
   "Wachter($t=r$)", and the free Jacobi (liberation) process (Demni–Hamdi–Hmidi) at matched $T_2$.
5. **Larger $N'$.** Check whether the $x$-independent excess keeps growing (§3, §4) and how the size distribution
   $w_k$ evolves. $N'=16$ needs about 8 GB for the BPS basis; the size decomposition is dense
   $D\times D$ ($D=12870$, about 2.6 GB per copy), so it needs a memory plan first.

## 9. Reproducibility

- Seeds: SYK 0–3 (parity family, size decomposition), pairs 100:101 … 106:107 (two-model); Haar 0–1, 200:201.
- $q=3$, $P=\lfloor N'/2\rfloor$, $S=\{N'-1,0,\ldots,s-2\}$.
- Code state: commit `081877d` + working tree (new files listed above).
- Memory: all runs under `scripts/memwatch.py`, peaks 0.01–0.86 GB (`memwatch.jsonl`).

## 10. Exact two-copy structure: transmission spectrum and the two-model fourth moment [derived + verified]

*Added 2026-09-29. Code: `src/size_decomposition.py` (`transmission_matrix`, `transmissions_from_weights`,
`johnson_idempotent_values`, `character_coeffs`); analysis `scripts/analyze_two_copy.py` → `two_copy_tables.md`;
test `test_transmissions_brute_force`. Extra data: 12 more independent pairs at $N'=10$ and $N'=14$ (seeds
110–133), size weights for $N'=14$ seeds 2–3 plus Haar. Peaks ≤ 0.7 GB under memwatch.*

### 10.1 The transmission spectrum
Define the superoperator $\Phi(X)=E[PXP]$ (average over the ensemble). By $U(N)$ invariance, $\Phi$ commutes with
$\mathrm{Ad}_g$ on $\mathrm{End}(\Lambda^P)$. The decomposition is multiplicity-free, so by Schur $\Phi$ acts on each
$V_k$ as a scalar $\phi_k$: **the transmission of a size-$k$ operator through the BPS projection.** Equivalently
$\phi_k=E\,\mathrm{Tr}[(P\otimes P^T)\Pi_k]/\dim V_k$, and for any operator $X$,
$E\,\mathrm{Tr}(PXPX^\dagger)=\sum_k\phi_k\|X^{(k)}\|^2$.

For two independent models,
$$E\,\mathrm{Tr}(P_1P_2P_1P_2)=E_1\,\mathrm{Tr}\big(P_1\,\Phi(P_1)\big)=\sum_k\phi_k\,E\|P_1^{(k)}\|^2
\;\Rightarrow\;\boxed{\,E\,\tau(P_1P_2P_1P_2)=\sum_k\phi_k\,w_k\,}$$

### 10.2 $\phi$ is a universal linear function of $w$ [derived]
**Realignment.** $P\otimes P^T$, the superoperator $X\mapsto PXP$, is the realignment
$\mathcal R\big(|P\rangle\rangle\langle\langle P|\big)$ of the rank-one operator whose twirl gives $w_k$. For
$Y=\sum_i|A_i\rangle\rangle\langle\langle B_i|$, $\mathcal R(Y)$ is the superoperator $X\mapsto\sum_iA_iXB_i^\dagger$.
$\mathcal R$ is $U(N)$-equivariant: $\mathrm{Ad}_g\,\mathcal R(Y)\,\mathrm{Ad}_g^{-1}=\mathcal R(\mathrm{Ad}_gY\mathrm{Ad}_g^{-1})$.
So it maps the commutant $\mathrm{span}\{\Pi_k\}$ to itself and commutes with twirling:
$$\phi_j=\sum_k\mathcal A_{jk}\,c_k,\qquad c_k=\frac{\|P^{(k)}\|^2}{\dim V_k}=\frac{d\,w_k}{\dim V_k},\qquad
\mathcal R(\Pi_k)=\sum_j\mathcal A_{jk}\Pi_j .$$
This holds **per realization** (for the $U(N)$-twirled $\phi$), not only in expectation. $\mathcal A$ depends only
on $(N,P)$.

**Computing $\mathcal A$ via the torus.**
1. $\mathrm{Tr}[\mathcal R(Y)\,\mathrm{Ad}_g]=\langle\langle\rho(g)^\dagger|Y|\rho(g)^\dagger\rangle\rangle$. So column $k$
   of $\mathcal A$ is the character expansion of the class function $g\mapsto\|\rho(g)^{(k)}\|^2$ (the size-$k$
   weight of the group element itself).
2. On the maximal torus, $\rho(g)=\mathrm{diag}(z^S)$ is diagonal. The size-$k$ part of a diagonal operator is its
   projection onto the $k$-th eigenspace of the Johnson scheme $J(N,P)$: the weight-zero subspace of $V_k$ carries
   the $S_N$ irrep $(N-k,k)$. So $\|\rho(g)^{(k)}\|^2=\sum_{S,T}(E_k)_{ST}z^S\bar z^T$, with $E_k$ the $k$-th
   Johnson idempotent, whose entries $e_k(r)$ depend only on $r=|S\setminus T|$.
3. With $m_r=\sum_{A\cap B=\emptyset,|A|=|B|=r}z^A\bar z^B$ one has
   $\sum_{|S\setminus T|=r}z^S\bar z^T=\binom{N-2r}{P-r}m_r$ and
   $\chi_{V_k}=|e_k|^2-|e_{k-1}|^2=\sum_r\big[\binom{N-2r}{k-r}-\binom{N-2r}{k-r-1}\big]m_r$.
   Matching coefficients of $m_r$ is a triangular system for each column.

**Verified:** brute force at $N=6$ (full $400\times400$ superoperator $P\otimes P^T$ projected with the Casimir
eigenbasis) agrees to $10^{-9}$ for two SYK and one Haar projector.

(Only the diagonal correlations $E[P_{SS}P_{TT}]$ enter the torus form of $\phi$ directly. Estimating them from
single realizations is noisy; for example, $\phi_0$ comes out 0.705 or 0.906 against the exact $a=0.771$ at $N'=8$.
The realignment route uses the full matrix through $w_k$ and is exact per realization.)

### 10.3 Exact sum rules [derived + verified per realization]
| identity | reason | max deviation ($N'=8$–14) |
|---|---|---|
| $\phi_0=a$ | $\Pi_0$ is the identity direction: $\mathrm{Tr}(P^2)/D$ | $5\times10^{-13}$ |
| $\sum_k\phi_k\dim V_k=d^2$ | trace of $P\otimes P^T$ $=|\mathrm{Tr}P|^2$ | $6\times10^{-15}$ (relative) |
| $r=T_2=\phi_1/a$ (half filling) | $U=1-2n_i$ has sizes 0, 1 only; $\|U^{(1)}\|^2=D$ | $2\times10^{-13}$ vs $\sum_kw_k\hat\chi_k(U)$ |

The third row gives **two independent exact routes to the decoder variance** (one via $w$ and characters, one via
the transmission of one-body operators), and they agree. The second row implies
$\sum_{k\ge1}(\phi_k-a^2)\dim V_k=-a(1-a)$. So the dimension-weighted mean transmission excess is essentially zero:
an excess at small sizes must be balanced by a tiny deficit spread over the huge large-size sectors.

### 10.4 Freeness = flat transmission [derived + num]
- **Haar:** $\phi_k=a^2$ for all $k\ge1$ up to $O(D^{-2})$ (e.g. $0.2767=a^2$ at $N'=12$). Hence
  $E\tau(P_1P_2P_1P_2)=2a^2-a^3$ (free), and $r=a$.
- **SYK:** mean over 4 realizations; seed-to-seed spread $\sim10^{-3}$ ($\phi$ is strongly self-averaging).

| $N'$ | $\delta\phi_k=\phi_k-a^2$ for $k=1,2,3,4,\ldots$ |
|---|---|
| 8 | +0.0351, −0.0031, −0.0036, +0.0047 |
| 10 | +0.0659, +0.0108, −0.0039, −0.0017, +0.0036 |
| 12 | +0.0839, +0.0245, +0.0024, −0.0027, −0.0006, +0.0024 |
| 13 | +0.0809, +0.0274, +0.0063, −0.0008, −0.0011, +0.0008 |
| 14 | +0.0902, +0.0333, +0.0092, −0.0000, −0.0017, −0.0001, +0.0016 |

Small operators are transmitted through the BPS projection **more** than freeness allows, increasingly so as $N'$
grows. $\delta\phi$ decays with size and oscillates slightly around 0 at intermediate $k$ (the compensation in
§10.3).

### 10.5 The two-model fourth moment: prediction vs measurement [verified]
Using $\phi_0=w_0=a$ and $\sum_kw_k=1$,
$$E\,\tau(P_1P_2P_1P_2)=(2a^2-a^3)+\sum_{k\ge1}\delta\phi_k\,w_k\qquad\text{(exact).}$$
The measured value uses independent pairs. The control-variate estimate regresses on $\nu-a$, whose expectation is
exactly 0 (§5).

| $N'$ | predicted (single-model data) | $\sum_{k\ge1}\delta\phi_kw_k$ | measured, plain (pairs) | measured, control variate | pred / free |
|---|---|---|---|---|---|
| 8 | 0.73054 | −0.00058 | – | – | 0.9992 |
| 10 | 0.56399 | +0.00313 | 0.56498 ± 0.00146 (16) | 0.56406 ± 0.00014 | 1.0056 |
| 12 | 0.41606 | +0.00827 | 0.41613 ± 0.00149 (4) | 0.41674 ± 0.00030 | 1.0203 |
| 13 | 0.29481 | +0.01053 | – | – | 1.0370 |
| 14 | 0.29661 | +0.01233 | 0.29645 ± 0.00074 (16) | 0.29664 ± 0.00014 | 1.0434 |
| Haar (all) | $=2a^2-a^3$ | 0.00000 | agrees | – | 1.0000 |

- Agreement is within $0.5\sigma$ at $N'=10$ and 14 (16 pairs each).
- At $N'=12$ (4 pairs) the plain mean agrees; the control-variate estimate is $2.3\sigma$ off, but with $n=4$ its
  error bar is unreliable.
- The earlier $3\sigma$ tension at $N'=14$ (§4, 4 pairs) was sampling noise.
- **The intrinsic non-freeness of two independent SYK BPS spaces is now an exact, parameter-free function of the
  single-model size spectrum:** $-0.08\%$ ($N'=8$), $+0.6\%$, $+2.0\%$, $+3.7\%$, $+4.3\%$ ($N'=10$–14).

### 10.6 Interpretation [interp]
- **One object controls both kinds of non-freeness:** the transmission excess $\delta\phi_k$, the extent to which a
  size-$k$ operator survives a sandwich between BPS projections beyond the free value $a^2$.
  - *Decoder variance* (the $x$-dependent effect of §3): the probe $U$ is one-body, so only
    $\delta\phi_1$ enters, $r-a=\delta\phi_1/a$.
  - *Intrinsic non-freeness* (the $x$-independent effect of §3–4): the second projector $P_1$ is itself spread over
    sizes with weights $w_k$, so the excess is $\sum_k\delta\phi_kw_k$. It is large because both $w_k$ and
    $\delta\phi_k$ sit at small $k$.
- **Freeness means flat transmission** ($\delta\phi\equiv0$); a Haar projector cannot tell operator sizes apart.
- **Chord comparison.** A size-$k$ operator is like matter with $\Delta_k=k/p$, so one expects
  $\phi_k/a\approx\langle q^{2\Delta_kn}\rangle$ (BLY 4.8).
  - $N'=12$: $\phi_1/a=0.686$ vs 0.690; for $k=2$ the excess over the size-0 value matches (0.047 vs 0.045).
  - For $k\ge3$ the exact excess falls to ~0 (and slightly below) faster than the chord curve.
  - The double-scaled formula (always $\ge a_\lambda$) has no analogue of the compensation sum rule of §10.3.
    Qualitatively right, quantitatively approximate beyond $k=2$.

### 10.7 What remains for the actual decoder $T_4$
The decoder's $T_4=\tau(PP'PP')$ has $P'=UPU$ *correlated* with $P$, so it does not factorize into
$E[P^{\otimes2}]\cdot E[P^{\otimes2}]$. An exact expectation-value route exists but needs the **four-copy**
structure:
- Averaging $U$ over its $U(N)$ orbit means $n_i\to n_v=\bar\psi(v)\psi(v)$ for a Haar-random orbital $v$.
- $E_v[v^{\otimes4}\bar v^{\otimes4}]$ is a sum over the 24 permutations in $S_4$. So $E\,T_4$ is a fixed
  combination of 24 index contractions $\sum\mathrm{Tr}(PE_{a_1b_1}PE_{a_2b_2}PE_{a_3b_3}PE_{a_4b_4})$, all
  computable from $P$ and the superoperator $M$.
- This is a finite computation, not a blocker. It does not yet give an interpretable spectral decomposition like
  §10.1. The natural target is a "two-body transmission" of the pair $U\otimes U$ through $P\otimes P$ that
  interpolates between the two-model value ($x=0$) and 1 ($x=1$). That is the next step.

