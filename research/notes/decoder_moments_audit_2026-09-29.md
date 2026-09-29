# Audit of `decoder_moments_2026-09-29.md`

*2026-09-29. Claim-by-claim truthfulness check of the running note before it is written up in
`research/pdfs/progress_2026-09-29.pdf`.*

Legend: ✅ verified as stated · ⚠️ true but needed qualification or rewording (corrected in the report) ·
❌ not supported (moved to open questions) · ➕ new result found during the audit.

## Verification runs performed for this audit
- Odd moments at even $N'$: $\max|T_1|,|T_3|\le 8\times10^{-16}$ per realization ($N'=10,12,14$).
- Odd-size weights at even $N'$: $\max w_{\rm odd}\le 2\times10^{-28}$ ($N'=8$–14), i.e. exactly zero.
- $N'=13$ sum rule with the correct odd-$N'$ conversion $T_2=r(1-c^2)+c^2$, $c=1/13$: measured 0.6193,
  predicted 0.6176 ($-0.28\%$).
- Parity-family sum rule with 4 seeds ($N'=12,13,14$, all $s$): agreement within $2.5\sigma$ everywhere.
  Largest tensions: $+2.4\sigma$ ($N'=14$, $s=7$) and $-2.0\sigma$ ($N'=13$, $s=3$).
- Closed form $\hat\chi_k(U_{\rm refl})=1-4k(N+1-k)/(N(N+1))$: deviation $\le10^{-16}$ for $N=7,\ldots,20$.
  Also $c_P=\sum_{ab}E_{ab}E_{ba}=P(N+1-P)$ on $\Lambda^P$.

## Findings

| § of note | claim | verdict | action |
|---|---|---|---|
| 0.1, 1 | decoder law = principal-angle law of $B$, $B'=UB$ | ⚠️ | Holds for $\lvert\mu\rvert$; the sign of $\mu$ is extra. At half filling the spectrum is symmetric ($T_{\rm odd}=0$ verified), so the two are equivalent there. |
| 1 | $T_{\rm odd}=0$ at half filling | ✅ | $\le 8\times10^{-16}$ per realization |
| 2 | "moving the $U$'s together leaves $\mathrm{tr}(U^m)$ × signs" | ⚠️ | Only for even $m$ ($U^m=1$). For odd $m$ one is left with $\mathrm{tr}(U\cdot\text{word}')$. |
| 2 | Z₂ rule: weight $x$ per chord separating an odd number of $U$'s | ⚠️ | Exact *per chord*. Treating chords as independent (factorization) is the double-scaled approximation; at finite $N$ different chords share site indices. |
| 3 | $x_s$ formula; $x=0$ exactly at $s=N'/2$ | ✅ | test `test_x_zero_at_half` |
| 3 | "Haar controls $1.000\pm0.001$ in every column" | ⚠️ | True for the ratio columns, not the KS column |
| 3 | "residual explained exactly in §6.3" | ⚠️ | Predicted $r/a=1.006$ vs measured $1.010\pm0.002$ ($N'=14$): consistent at $2\sigma$. Not "exactly". |
| 3 | "94 % of the $T_4$ excess disappears as $x\to0$" | ⚠️ | The comparison is between two *different probes* (the decoder and the $s=N'/2$ parity), not a continuous deformation. Reworded. |
| 3 | $T_4,T_6$ exceed free at $x\approx0$, growing with $N'$ | ✅ | |
| 3 | $P_n$ extraction from $T_2(x_s)$ ill-posed | ✅ | |
| 4, 10.5 | two independent SYK BPS spaces are non-free (+0.6, 2.0, 4.3 %) | ✅ | 16/4/16 pairs; exact prediction agrees |
| 5 | $E[P]=a\mathbb 1$, $E[\nu]=a$ | ✅ | $U(N)$ invariance of the complex-Gaussian coupling law (checked in `q_scan.rand_form`) |
| 6.1 | $V_k$ structure, Casimir $2k(N+1-k)$, Krawtchouk characters | ✅ | tests |
| 6.2 | sum rule $E[T_2]=\sum_k\bar w_k\hat\chi_k$; $w_0=a$ exactly | ✅ | |
| 6.3 | table row $N'=13$ compares $\sum w\hat\chi$ with $r$ | ⚠️ | At odd $N'$, $r\ne T_2$. Corrected: 0.6176 vs 0.6193 ($-0.28\%$). Conclusion unchanged. |
| 6.3 | parity family "follows the same rule for every $s$" | ⚠️ | Within $2.5\sigma$ with 4 seeds; reworded. |
| 6.4 | only even sizes at half filling | ✅ numerically; derivation is a sketch | kept as [derived (sketch) + verified] |
| 6.5 | $w_{2n}\approx P_n$, $\hat\chi_{2n}\approx x^{2n}$ | ✅ as an empirical match (10–20 %) | $n\leftrightarrow k/2$ not derived |
| 6.5 | "termwise errors cancel, *which is why* BLY is good to 2 %" | ⚠️ | Causal claim not shown. Reworded to "consistent with". |
| 7, 0.8 | "free compression = zero-length truncation" falsified "in strong form" by SYK at $x\to0$ | ⚠️ | Reframed with a sharper, exact statement: already for **Haar** $P$, $T_4^{\rm free}=a^2+Q_4$ with $Q_4=a^2(1-a)$ coming entirely from nonzero sizes. So "free = zero-size truncation" fails identically beyond $m_2$, independent of SYK. The SYK statement (removing $x$ does not recover freeness) stands separately. |
| 10.2 | $\phi=\mathcal A\,c$ universal, per realization | ✅ | brute force $N=6$, $10^{-9}$ |
| 10.3 | "two independent exact routes to the decoder variance" | ⚠️ | Both give the $U(N)$-*orbit-averaged* variance of a realization, not the single-site value. Clarified. |
| 10.4 | Haar $\phi_k=a^2+O(D^{-2})$ | ✅ | Weingarten: $\phi_{k\ge1}=a(aD^2-1)/(D^2-1)$ |
| 10.6 | "one object ($\delta\phi$) controls both kinds of non-freeness" | ✅ as stated for $m_2$ and the two-model $m_4$ | Not a statement about $R$ |
| 11.1 | "$x=0$ corresponds to a Haar-random probe" | ⚠️ | Exact content: $E_g\,T_4(g)=\sum_k\phi_kw_k$ (per realization, $g$ Haar in $U(N)$) and $E_gT_2=a$. The correspondence with independent models holds in expectation. Labelled [interp]. |
| 11.2 | exact decomposition $T_4=a^2+2a(T_2-a)+S+R$ | ✅ | test |
| 11.3 | "about 79 % fixed exactly" | ⚠️ | 77–79 % at $N'=13,14$; 94 % at $N'=12$. Stated per $N'$. |
| 11.4 | size-diagonal and coherent-mixture closures fail | ✅ | |
| 11.4 | "$\Gamma$ concentrates on small two-copy irreps" | ❌ | Inferred to explain the failure, not measured. Moved to open questions as a conjecture. |
| 11.4 | $R/R(\mathbb 1)\sim\eta$ | ✅ as a rough observation, explicitly not a law | |

## ➕ New during the audit (derived + verified)
For the single-site reflection $U_v=1-2n_v$, averaging over the orbit ($v$ a Haar-random orbital) gives
$E_v[U_vXU_v]=X-\frac{2}{N(N+1)}C(X)$, with $C$ the adjoint Casimir. Hence
$$\hat\chi_k(U_{\rm refl})=1-\frac{4k(N+1-k)}{N(N+1)},\qquad
r\ (\text{half filling})=E[T_2]=1-\frac{4\,\langle k(N+1-k)\rangle_w}{N(N+1)} .$$
**The decoder variance measures the mean adjoint Casimir ("operator size") of the BPS projector.** Numerically
this matches $\sum_kw_k\hat\chi_k$ to all digits. The mean size $\langle k\rangle_w$ grows as 0.49, 0.83, 1.23,
1.64 for $N'=8,10,12,14$.

## Further corrections made while writing the report
Found by cross-checking every table entry in the report against the saved analysis tables
(`results/data/parity_family_2026-09-29/*_tables.md`, `summary.csv`, the size/two-model/T₄ jsonl files):
- **Two-model fourth moment at $N'=12$** (§4/§10.5 of the note): only 4 pairs. The plain mean
  (0.41613 ± 0.00149) agrees with the prediction 0.41606. The control-variate mean (0.41674 ± 0.00030) is
  0.16 % (≈2σ) above it. The report shows both. "Agrees to ≤0.02 %" holds only at $N'=10,14$ (16 pairs each).
- **$R/R(\mathbb 1)\sim\eta$** (§11.4): the data give 0.23–0.25 at $N'=13,14$, between $\eta^2\approx0.11$–0.14 and
  $\eta\approx0.33$–0.38. The report says "between", not "like η".
- **Casimir superoperator test** is at $N=6$ only (not 6, 7). The size-component reconstruction residual is
  tested at $<10^{-10}$ ($N=8$).
- **$a_\lambda$ vs $a$:** 2–8 % at $N'=8$–14, from the $P_0$ column of the size tables.
- **Notation in the report:** the fermion-number sector is called $F$ (the code and notes use $P$), so that $P$
  always means the BPS projector.
