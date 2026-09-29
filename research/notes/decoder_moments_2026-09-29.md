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
   §5), but $E[\nu^2]$ exceeds the free value by $1.4, 2.0, 3.4\%$ at $N'=10, 12, 14$. This is the $x$-independent part
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
6. **[falsified, strong form] "Wachter = zero-length truncation of the chord computation for all moments"**
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
| 10 | 1.009 ± 0.006 | 1.014 ± 0.005 | 1.017 ± 0.004 | 0.026 / 0.010 |
| 12 | 0.999 ± 0.003 | 1.021 ± 0.004 | 1.040 ± 0.004 | 0.032 / 0.007 |
| 14 | 0.993 ± 0.002 | 1.034 ± 0.003 | 1.077 ± 0.004 | 0.038 / 0.007 |

- **Independent SYK BPS spaces are not free relative to each other,** and the deviation grows with $N'$.
- **They reproduce the $x=0$ parity numbers at $N'=10,12$.** At $N'=14$ the parity $x=0$ decoder is somewhat
  more non-free ($T_4$: 1.056 vs 1.034). Its residual correlation is the nonzero $\hat\chi_{k\ge2}(U_{N/2})$ of §6.3.
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
| 14 (2 seeds) | 0.4248 | 0.3707 | 0.1629 | 0.0416 | 0 | 0.6371 | 0.6420 ± 0.0037 |

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

1. **Exact $T_4$ from two-copy size data [derived route, to do].** For two independent models,
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
