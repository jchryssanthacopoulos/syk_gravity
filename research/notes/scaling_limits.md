# Scaling limits of the SYK uplift decoder: where the Wachter law lives

*Short companion note to* Fragile Fortuity: The Gravity of BPS Overlaps in Supersymmetric SYK. *Distilled from the working notes §§2–3, 6–7; every $a,b,\lambda_\pm$ below is checked against the exact BPS-dimension formula (eq. 125) and exact diagonalization of the one-flavor $\mathcal N{=}2$ wedge model.*

## Setup: the two Wachter parameters and the band edges

The uplift decoder $\mathcal O_T=P_B\,\Pi_T\,P_B$ has eigenvalues $\lambda=\cos^2\theta$, the squared cosines of the principal angles between the BPS space $\mathcal B$ (dimension $d$) and the slot $V$ (dimension $m$) inside the enlarged Fock space (dimension $D$). To leading order in $1/D^2$ this is a **Wachter (Jacobi / MANOVA) ensemble**, controlled by exactly two aspect ratios:

$$a=\frac dD\quad(\text{BPS fraction}),\qquad b=\frac mD\quad(\text{slot fraction}).$$

The Wachter law has a bulk band with two hard edges plus possible atoms at the endpoints $\lambda=0,1$:

$$\lambda_\pm=\Big(\sqrt{a(1-b)}\pm\sqrt{b(1-a)}\Big)^2,\qquad
p_0=\max\!\Big(0,1-\tfrac{b}{a},\,1-\tfrac{1-b}{1-a}\Big),\ \dots$$

Everything that follows is what happens to $a$, $b$, and hence $\lambda_\pm$, in the three scaling limits. The single structural fact that makes the story clean is that the model sits at **half filling**, where $b$ is pinned.

## $b$ is elementary and exact — it is $\tfrac12$ in every limit

$b$ needs no scaling argument. With the enlarged charge sector $M=N+1$ and occupation $p$,

$$b=\frac mD=\frac{\binom Np}{\binom{N+1}{p}}
=\frac{(N+1-p)!/(N-p)!}{(N+1)!/N!}=1-\frac{p}{N+1}
\ \xrightarrow{\,p=(N+1)/2\,}\ \tfrac12 .$$

This is exact at every $N,p,q$; at half filling $b=\tfrac12$ identically. So the whole family lives on the line $b=\tfrac12$, and the edge formula collapses to a one-parameter expression in $a$ alone:

$$\boxed{\,b=\tfrac12\ \Rightarrow\ u\equiv a+b-2ab=\tfrac12\ (\text{independent of }a),\quad
\lambda_\pm=\tfrac12\pm\sqrt{a(1-a)},\quad
p_0=p_1=\max\!\big(0,\,1-\tfrac1{2a}\big).\,}$$

*(Check: $\big(\sqrt{a\cdot\tfrac12}\pm\sqrt{\tfrac12(1-a)}\big)^2=\tfrac12\big(1\pm2\sqrt{a(1-a)}\big)$.)*

Atoms appear only when $a>b=\tfrac12$. All the dynamics is therefore carried by the single number $a=d/D$, and the three limits are three behaviors of $a$.

## $a$ from the BPS count (eq. 125)

The BPS states are $\pm1$ lattice walks confined to the fortuity strip $|2p-N|\le q-1$; the exact count is the method-of-images alternating sum

$$d=\dim\mathcal B_o^p=\widetilde D_p-\widetilde D_{N-p-q},\qquad
\widetilde D_p=\sum_{n\ge0}(-1)^n\binom{N}{p-nq}.$$

Normalizing by the leading (central, $n=0$) binomial and using the local CLT
$\binom{N}{p-k}\big/\binom Np\to e^{-2k^2/N}$ near half filling (with $k=nq$), the alternating image sum telescopes into a Jacobi theta:

$$\boxed{\,a=\frac dD\ \longrightarrow\ \sum_{n=-\infty}^{\infty}(-1)^n e^{-2q^2 n^2/N}
=\vartheta_4\!\big(0,\mathfrak q\big),\qquad \mathfrak q=e^{-2q^2/N}\ (\text{the DSSYK nome}).\,}$$

The behavior of this theta in the three limits is the whole story. The controlling ratio is $\eta\equiv q^2/N=(\text{strip width}/\text{diffusion length})^2$.

---

## The three limits, and where Wachter lives

| limit | condition | $a=d/D$ | $b$ | band | verdict |
|---|---|---|---|---|---|
| **single-scaled** | $q$ fixed, $N\to\infty$ | $\to 0$ (rate $r(q)$) | $\tfrac12$ | $\to\delta(\lambda-\tfrac12)$ | **not Wachter** (band vacuous) |
| **double-scaled** | $\eta=q^2/N$ fixed | $\vartheta_4(0,e^{-2\eta})=O(1)$ | $\tfrac12$ | genuine $O(1)$ two-edge band | **pure Wachter** ($+\,1/D^2$ handles) |
| **triple-scaled** | double-scaling $+$ zoom onto near-BPS edge | $\lambda=E/(E+\sigma_a^2)$ | — | resolves near-extremal throat | **the deviation from Wachter** |

### 1. Single-scaled ($q$ fixed, $N\to\infty$): the band is vacuous

Here $\eta=q^2/N\to0$, so $\mathfrak q=e^{-2\eta}\to1$ and the theta degenerates. Directly from the count, confinement to a fixed-width strip is exponentially rare:

$$a\sim e^{-r(q)N},\qquad r(q)=-\ln\cos\tfrac{\pi}{2q}\approx\frac{\pi^2}{8q^2}\ (\text{large }q),$$
$$q=3:\ r=\ln2-\tfrac12\ln3=0.1438.$$

With $a\to0$ and $b=\tfrac12$ fixed:

$$\lambda_\pm=\tfrac12\pm\sqrt{a(1-a)}\ \longrightarrow\ \tfrac12\pm\sqrt a\ \longrightarrow\ \tfrac12,
\qquad \text{band width }=2\sqrt{a(1-a)}\to0,$$

so the bulk collapses to $\delta(\lambda-\tfrac12)$. Since $a<b$, there are **no atoms** ($p_0=p_1=0$): fortuity's protected core is a vanishing finite-$N$ resonance. The Wachter law is *formally present but content-free* — its band is a single point, and the entire nontrivial spectrum is the near-BPS tail that sits outside it. **This is not a Wachter limit** in any useful sense; the band is a delta and the tail is everything.

### 2. Double-scaled ($\eta=q^2/N$ fixed): pure Wachter

Now $\mathfrak q=e^{-2\eta}$ is $O(1)$, so $a=\vartheta_4(0,e^{-2\eta})$ is a genuine $O(1)$ number strictly between $0$ and $1$. Both aspect ratios are $O(1)$:

$$a=\vartheta_4(0,e^{-2\eta})\in(0,1),\qquad b=\tfrac12,$$
$$\lambda_\pm=\tfrac12\pm\sqrt{a(1-a)}\quad(\text{two hard edges, both away from }0,1\text{ when }a\ne\tfrac12),$$
$$p_0=p_1=\max\!\big(0,1-\tfrac1{2a}\big).$$

This is the home of **plain Wachter**: the overlap ("decoder") chords have crossing weight $1/D^2=e^{-2S_0}$, which vanishes automatically at large $N$, so the ensemble is at the *free point* — non-crossing dominates and one gets the plain (not $q$-deformed) Wachter disk, dressed only by a genus / $1/D^2$ handle expansion with $O(1)$ parameters. The DSSYK nome $\mathfrak q$ enters *only* through the parameter $a=\vartheta_4(\mathfrak q)$, not through a $q$-deformation of the law.

The atoms turn on at a sharp transition: $p_0>0\iff a>b=\tfrac12$, i.e.

$$\vartheta_4(0,e^{-2\eta^\star})=\tfrac12\ \Rightarrow\ \boxed{\eta^\star=0.685}.$$

For $\eta>\eta^\star$ the protected sector ($\lambda=0,1$ atoms of weight $1-\tfrac1{2a}$) **survives as $N\to\infty$** — fortuity's protected directions stabilize into an $O(1)$ fraction, dimensionally forced by $d>m$. $\eta^\star$ is simultaneously an **Airy$\leftrightarrow$Bessel edge-universality transition**: soft (Airy) edge at $\lambda_-$ with no atoms below, hard (Bessel) edge with atoms above.

**Caveat (the one open clause).** "Pure Wachter" here is the band/disk law, which is provable (both parameters $O(1)$, Weingarten + eq. 125 give a Wachter+handles series to all orders in $1/D^2$). *Exact* Wachter additionally needs the freeness order parameter $t_2=m_2^{\rm SYK}-m_2^{\rm W}\to0$ — the near-BPS continuum must stay **starved** so the soft-edge tail does not leak. That piece is non-perturbative in $1/D$ (invisible to the genus series) and is exactly item 3.

### 3. Triple-scaled (double-scaling $+$ near-BPS zoom): the part that deviates

Take the double-scaled model and zoom onto the near-BPS edge $\lambda\to0$, resolving the near-extremal throat — the individual near-BPS states, the gap, $\rho_H$. The exact relation between a decoder eigenvalue and the underlying near-BPS Rayleigh quotient $E$ (verified to $10^{-16}$) is

$$\lambda=\frac{E}{E+\sigma_a^2},\qquad
\sigma_a^2\sim\binom{N/2}{q-1}\sim\frac{N^2}{8}\ (\text{the zoom factor}),$$

so the tail sits at $\lambda\sim1/\sigma_a^2$. Changing variables,

$$\rho_{\rm tail}(\lambda)=\frac{\sigma_a^2}{(1-\lambda)^2}\,\rho_H\!\Big(\frac{\sigma_a^2\lambda}{1-\lambda}\Big),$$

and the **triple-scaling limit** holds $\eta$ fixed while zooming $\lambda\to0$ with $E\equiv\sigma_a^2\lambda/(1-\lambda)$ fixed, giving exactly

$$\rho_{\rm tail}(\lambda)\,d\lambda\ \longrightarrow\ \rho_H(E)\,dE .$$

Three windows: the band is $\lambda=O(1)\Leftrightarrow E\sim\sigma_a^2$; the tail is $\lambda=O(1/\sigma_a^2)\Leftrightarrow E=O(1)$. In modular language the band is the theta *winding* (chord) rep and the tail is its *momentum* (energy) dual — triple-scaling is the passage to the dual frame.

The point for Wachter: **triple-scaling zooms into precisely the part of the spectrum that deviates from the Wachter law.** The band (item 2) is where the law holds; the near-BPS edge it magnifies is where it fails — the macroscopic soft-edge leak below $\lambda_-$ ($\sim9\%$ weight, power law $\rho\sim\lambda^{-0.6}$ vs a sharp Airy edge for generic Haar), which is the small-$E$ (near-kernel) structure of $\rho_H$ imaged into overlaps.

> **Correction carried from §7 of the notes.** $E$ is a *Rayleigh quotient* (harmonic-admixture measure), **not** an energy eigenvalue: numerically the tail states have $E\sim10^{-2}$, far below the smallest nonzero $H_N$ eigenvalue ($\approx2.25$), with occupied component $\gtrsim99.99\%$ old-BPS harmonic. So $\rho_H$ is the admixture distribution — the small-$\sigma_b^2$ near-kernel of the compressed contraction $\iota_{\bar D}^\dagger\iota_{\bar D}\big|_{\mathcal B_N^{p-1}}$ (itself a MANOVA object) — **not** an energy-continuum DOS leaking through a gap. Whether its non-generic (beyond-Haar) excess is $O(1)$ (genuine tail, fortuity survives) or $O(1/D^2)$ (pure Wachter) in the DS limit is the one unresolved free-probability question.

---

## One-line summary

$b=\tfrac12$ always; the physics is entirely in $a=\vartheta_4(0,e^{-2q^2/N})$ and $\lambda_\pm=\tfrac12\pm\sqrt{a(1-a)}$. Single-scaling sends $a\to0$ and the band to a point (Wachter vacuous); double-scaling holds $a=O(1)$ and delivers the genuine two-edge Wachter law with protected atoms above $\eta^\star=0.685$; triple-scaling magnifies the near-BPS edge $\lambda=E/(E+\sigma_a^2)$ — which is exactly where the Wachter law breaks.
