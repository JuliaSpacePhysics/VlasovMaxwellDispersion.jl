# Derivation

The susceptibility of a hot magnetized plasma is **one helical-orbit phase-space integral**. 
This document shows how this integral collapses, step by step, under specific assumptions.

## Conventions

Normalize `v → v/c`, `ω → ω/Ω_ref`, `k → k c/Ω_ref`.

Let $Ω=qB/(γmc)$ be the relativistic gyrofrequency, $Ω₀=γΩ$ the signed constant rest
gyrofrequency, and $Π ≡ ω_p/Ω_{ref}$.

Phase-space measure: `d³p = 2π p⊥ dp⊥ dp∥ = 2π γ dγ dp∥`.

## 1. Linear Vlasov → orbit integral

Linearize relativistic Vlasov–Maxwell about a gyrotropic `f₀(p⊥,p∥)`. 
With `θ=0` the orbit phase is

    β = −z[sin(φ + Ω τ) − sin φ] + (ω − k∥ v∥)τ.

Set `s=Ωτ`, and define the three dimensionless orbit parameters

    a ≡ (ω − k∥ v∥)/Ω = (γω − k∥p∥)/Ω₀,   z ≡ k⊥ v⊥/Ω = k⊥ p⊥/Ω₀,

with `φ` the azimuthal angle of `p`. Every field component then reduces to the single
orbit integral

    g(φ,z) = Ω⁻¹ ∫₀^∞ exp[−i z sin(φ+s) + i a s] ds.

## 2. Gyrophase average → susceptibility tensor

Assembling the current moment from `g(φ,z)` by a gyrophase average gives

    χ = (Π²)/(ω Ω₀) ∫ d³p · S

with

$$
S = -i\Omega
\begin{pmatrix}
p_⊥ U G_{33} & p_⊥ U G_{32} &
  p_⊥ ∂_{p_∥}f_0 G_{31} - p_⊥ V G_{33} \\
p_⊥ U G_{23} & p_⊥ U G_{22} &
  p_⊥ ∂_{p_∥}f_0 G_{21} - p_⊥ V G_{23} \\
p_∥ U G_{13} & p_∥ U G_{12} &
  p_∥ ∂_{p_∥}f_0 G_{11} - p_∥ V G_{13}
\end{pmatrix},
$$

the operators

$$
U = ∂_⊥f₀ + \frac{k_∥}{ωγ}(p_⊥ ∂_∥f₀ − p_∥ ∂_⊥f₀),
\qquad
V = \frac{k_⊥}{ωγ}(p_⊥ ∂_∥f₀ − p_∥ ∂_⊥f₀),
\\

$$

and the gyrophase matrix `G[g]`

$$
G_{i j} \equiv \frac{1}{2 \pi} \int_0^{2 \pi} d φ\, e^{i z \sin φ }\left(\begin{array}{ccc}
g & i \frac{∂ g}{∂ z} & \frac{i}{z} \frac{∂ g}{∂ φ} \\
g \sin φ & i \frac{∂ g}{∂ z} \sin φ & \frac{i}{z} \frac{∂ g}{∂ φ} \sin φ \\
g \cos φ & i \frac{∂ g}{∂ z} \cos φ & \frac{i}{z} \frac{∂ g}{∂ φ} \cos φ
\end{array}\right).
$$

`S` has identical form in physical or normalized variables; the per-species `m_s, q_s`
sit only in the `Ω₀, Π` prefactor.

The `V` dependence (last column of `S`) factors into a compact closed form, the `e∥e∥`-Bernstein term (§3).

From here the two evaluators diverge: **A** uses the periodicity of `g`; **B** expands
it in harmonics.

## 3. Closing the orbit integral

**(A) Finite-period reduction (Qin's complex-order form)** Gyrophase symmetry `g(φ,z)=g(φ+2π,z)` splits the semi-infinite `s` integral into gyroperiods,

    g(φ,z) = c₀/(2π) ∫₀^{2π} exp[-i z sin(φ+η) + i a η] dη,
    c₀ = −e^{-iπa} π / (i Ω sin πa),

retaining every resonance in the `1/sin(πa)` factor (zeros at `ωγ−k∥p∥=nΩ₀`). After some algebra, we can reduce `G₁₁ = −π J_{−a}(z)J_a(z)/(iΩ sin πa)`, and with

    σ0 ≡ π J_{−a}(z)J_a(z)/sin πa,   Q ≡ a σ0,   Q' = a ∂_z σ0,

the resonant 3×3 `T(a,z)` closes in `Q, Q', J'_{−a}J'_a, a, z` and `r=p∥/p⊥`:

    T = [
      a/z² (Q−1)              i/(2z) Q'                 (Q−1) r/z
     −i/(2z) Q'               π J'_{−a}J'_a/sin πa + a/z²   −i r Q'/(2a)
      (Q−1) r/z               i r Q'/(2a)                (Q/a) r²
    ].

(Qin Eq. 35's printed `T₁₂,T₂₁,T₂₃,T₃₂` have typos.) Folding `V` in,

    S = e∥e∥ (Ω/ω)(p∥/p⊥)(p⊥ ∂∥f₀ − p∥ ∂⊥f₀) + p⊥ U T.

The susceptibility splits into a resonant and a non-resonant (Bernstein) term:

$$
\chi = \chi_{T} + \chi_{B}, \qquad
\chi_{T} ≡ \frac{\Pi^2}{\omega \Omega_0} \int d^3p\; p_\perp U T(a,z),\quad
\chi_B ≡ \frac{\Pi^2}{\omega \Omega_0} \int d^3p\;\frac{\Omega}{\omega}\left(\frac{1}{p_\|}\frac{\partial f_0}{\partial p_\|} - \frac{1}{p_\perp}\frac{\partial f_0}{\partial p_\perp}\right)p_\|^2 .
$$

Price: `a` is complex ⇒ complex-order Bessel; and relativistically `sin πa=0` traces
curves in the momentum plane, leaving a 2-D integral plus residues (§3.1) — versus the
clean 1-D reduction B reaches in the separable case.

See [A_evaluator.md](A_evaluator.md) for more details.

**(B) Harmonic sum.** Insert `e^{iz sinφ}=Σ_m Jₘ(z)e^{imφ}`; the orbit integral becomes
harmonic sums like `Σₙ Jₙ(z)²/(a−n)`, and

$$
S = \mathbf{e}_{\|} \mathbf{e}_{\|}\, \frac{Ω}{ω}\left(\frac{1}{p_{\|}} \frac{∂ f_0}{∂ p_{\|}}-\frac{1}{p_{\perp}} \frac{∂ f_0}{∂ p_{\perp}}\right) p_{\|}^2  +\sum_{n=-\infty}^{\infty} \frac{Ω\, p_{\perp} U}{ω-k_{\|} v_{\|}-n Ω}\, \mathbf{T}_n,
$$

with the per-harmonic tensor

$$
\mathbf{T}_n \equiv\left(\begin{array}{ccc}
\frac{n^2 J_n^2}{z^2} & \frac{i n J_n J_n^{\prime}}{z} & \frac{n J_n^2 p_{\|}}{z p_{\perp}} \\
-\frac{i n J_n J_n^{\prime}}{z} & \left(J_n^{\prime}\right)^2 & -\frac{i J_n J_n^{\prime} p_{\|}}{p_{\perp}} \\
\frac{n J_n^2 p_{\|}}{z p_{\perp}} & \frac{i J_n J_n^{\prime} p_{\|}}{p_{\perp}} & \frac{J_n^2 p_{\|}^2}{p_{\perp}^2}
\end{array}\right),
\qquad
𝓣_n \equiv p_⊥^2 \mathbf{T}_n .
$$

Using the ring kernel $R_n \equiv (n/z)J_n=\tfrac12(J_{n-1}+J_{n+1})$, $𝓣_n \equiv p_\perp^2 T_n$ reads

$$
𝓣_n=
\begin{pmatrix}
p_\perp^2 R_n^2 & i\,p_\perp^2 R_n J_n' & p_\parallel p_\perp R_n J_n\\[2pt]
-i\,p_\perp^2 R_n J_n' & p_\perp^2 J_n'^2 & -i\,p_\parallel p_\perp J_nJ_n'\\[2pt]
p_\parallel p_\perp R_n J_n & i\,p_\parallel p_\perp J_nJ_n' & p_\parallel^2 J_n^2
\end{pmatrix}.
$$

Note: For *real* $(p_\|,p_\perp)$, $𝓣_n$ is the Hermitian
outer product $\mathbf{v}_n\mathbf{v}_n^\dagger$ of $\mathbf{v}_n=(p_\perp R_n,-i p_\perp J_n',p_\|J_n)$.

The susceptibility splits into harmonic contributions plus one non-resonant term:

$$
\chi = \sum_n \chi_n + \mathbf{e}_\| \mathbf{e}_\|\, \chi_B,\qquad
\chi_n = \frac{\Pi^2}{\omega \Omega_0} X_n,\quad
X_n ≡ \int d^3p\;\frac{p_\perp U \Omega}{\omega - k_\parallel v_\parallel - n\Omega}\;\mathbf{T}_n = \int d^3p\;\frac{p_\perp U}{a-n} \mathbf{T}_n.
$$


Keep `|n| ≤ nmax ≈ k⊥ρ`; the sum converges slowly for large `z`. The payoff: each raw
`X_n` **factorizes when `f₀` is separable** (§4), which is exactly what produces the two
1-D primitives. 

(`σ0 = Σ_n Jₙ(z)²/(a−n) = π J_{−a}J_a/sin πa` is the Lerche–Newberger Bessel identity linking A and B.)

## 4. The momentum integral: where the primitives appear

Either closure leaves the same shape,

    χ = (Π²/ω Ω₀) ∫ d³p [ e∥e∥-Bernstein + (p⊥ U) · K ],

with kernel `K` = closed `T(a,z)` (A) or harmonic sum `Σ_n T_n/(a−n)` (B). 


With 𝒰 ≡ (ωγ/p⊥) U = k∥ ∂∥f₀ + (ωγ − k∥p∥)/p⊥ · ∂⊥f₀ = ω ∂f₀/∂γ + k∥ ∂f₀/∂p∥ , the harmonic contribution is,

$$
X_n = 2π ∫ dp_⊥ dp_∥ \; U 𝓣_n /(a − n) = (2π/ω) ∫ dγ dp_∥ \;𝒰 𝓣_n /(a - n)
$$

**(0) Relativistic, coupled `f₀(p⊥,p∥)`.** Work in `(γ,p∥)`. The denominator is nonlinear
in `p∥` through `γ`, but at fixed `γ` it linearizes:

    ω − k∥v∥ − nΩ = (ωγ − k∥p∥ − nΩ₀)/γ = −(k∥/γ)(p∥ − ζ_n(γ)),   ζ_n(γ) = (ωγ − nΩ₀)/k∥,

a clean rational pole in `p∥`:

    χ_n = (Π²/ω²)·(−2π/k∥) ∫ dγ ∫_{|p∥|<√(γ²−1)} dp∥ · 𝒰 𝓣_n /(p∥ − ζ_n(γ)).

This straightens A's `sin πa=0` resonance curve into the line `p∥=ζ_n(γ)`. What stays
coupled is only `z=(k⊥/Ω₀)√(γ²−1−p∥²)`, so the `p∥` integral is the analytic `𝒞` branch
over a finite interval — no clean primitive, **still 2-D**.

**(3) Maxwellian / Maxwell–Jüttner / cold.** Gaussian `f` closes both primitives in
*closed form*: `𝒞→Z` (plasma dispersion function), `P⊥→Γ_n=Iₙ(λ)e^{−λ}`. No quadrature
— fastest. (Still evaluator B: the `n`-sum is truncated; only per-`n` moments are closed.)

## 5. The outer-coordinate density `I`

The cyclotron resonance is a **pole in $p_\parallel$** at $p_\parallel=\zeta_n$, at fixed outer
coordinate. Integrating inner $p_\parallel$ first — continued past the pole (Plemelj
boundary value + growing-sheet Landau residue, §5.3) — leaves a regular outer-coordinate density $I$:

$$\chi=\frac{\Pi^2}{\omega^2}\int I(\xi)\,d\xi.$$

Closures A and B differ only in how they reach the poles: B sums one single-pole transform per
harmonic (truncating the tail), A peels the in-range poles from the closed kernel and sums the rest.

The following reduction is **not** closure-general: the non-relativistic $I(p_\perp)$ factors into parallel moments
× perp Bessel weights (§5.1) only because B's $T_n$ cleanly separates $p_\parallel$-polynomials from the
$z(p_\perp)$-Bessel factors. A's complex-order kernel ties $a(p_\parallel)$ to $z(p_\perp)$, so it never
factors.

### 5.1 Non-relativistic

In $(p_\perp,p_\parallel)$ coordinates

$$
I = I(p_\perp)=\sum_n I_n(p_\perp)+\mathbf e_\parallel\mathbf e_\parallel\,I_B(p_\perp),
$$
$$
I_n(p_\perp)=2\pi\!\int dp_\parallel\;
\frac{p_\perp\mathcal U\; 𝓣_n}{\omega-k_\parallel p_\parallel-n\Omega},
\qquad
I_B(p_\perp)=2\pi\!\int dp_\parallel\;\Bigl(p_\perp p_\parallel\,\partial_\parallel f_0-p_\parallel^2\,\partial_\perp f_0\Bigr).
$$


With $\gamma=1$, three things happen: the pole $\zeta_n=(\omega-n\Omega)/k_\parallel$ loses its $\gamma$-dependence,
the Bessel argument $z=k_\perp p_\perp/\Omega$ depends on $p_\perp$ alone, and $𝓣_n$
is polynomial in $p_\parallel$. Now $p_\perp$ and $p_\parallel$ decouple, so

Unlike $I(\gamma)$, $I(p_\perp)$ **does** factor into parallel moments and perp Bessel weights.

**Parallel Landau moments.** 
With the recurring Landau combination ($p_\perp\mathcal U=(\omega-k_\parallel p_\parallel)\,\partial_\perp f_0 +k_\parallel p_\perp\,\partial_\parallel f_0$), the moment of the full numerator $p_\perp\mathcal U$ reduces to a linear combination of the two slices:

$$
D_m(p_\perp)\equiv\int\!\frac{p_\parallel^{\,m}\,p_\perp\mathcal U}{\omega-k_\parallel p_\parallel-n\Omega}\,dp_\parallel
=\omega\,M^m_F-k_\parallel M^{m+1}_F+k_\parallel p_\perp\,M^m_T, 
$$

where

$$
M^m_F(p_\perp)\equiv\int\!\frac{p_\parallel^{\,m}\,\partial_\perp f_0}{\omega-k_\parallel p_\parallel-n\Omega}\,dp_\parallel,
\qquad
M^m_T(p_\perp)\equiv\int\!\frac{p_\parallel^{\,m}\,\partial_\parallel f_0}{\omega-k_\parallel p_\parallel-n\Omega}\,dp_\parallel,
$$

$$
M^m_F=-\frac1{k_\|}\mathcal C[p_\|^m\,\partial_\perp f_0](\zeta_n),\qquad
M^m_T=-\frac1{k_\|}\mathcal C[p_\|^m\,\partial_\| f_0](\zeta_n).
$$

Note without separability, the moments are recomputed at every $p_\perp$, otherwise it would pull $f_\perp(p_\perp)$ out of $M^m$ and kill the outer integral.

**Bessel weights.** The Bessel factors leave the $p_\parallel$ integral untouched; each power of $p_\parallel$ raises the
moment index ($p_\parallel^0\!\to\!D_0$, $p_\parallel^1\!\to\!D_1$, $p_\parallel^2\!\to\!D_2$).

Assembled, the resonant harmonic block reads

$$
I_n(p_\perp)=2\pi
\begin{pmatrix}
p_\perp^2 R_n^2\,D_0 & i\,p_\perp^2 R_n J_n'\,D_0 & p_\perp R_n J_n\,D_1\\[2pt]
-i\,p_\perp^2 R_n J_n'\,D_0 & p_\perp^2 J_n'^2\,D_0 & -i\,p_\perp J_nJ_n'\,D_1\\[2pt]
p_\perp R_n J_n\,D_1 & i\,p_\perp J_nJ_n'\,D_1 & J_n^2\,D_2
\end{pmatrix},
$$

**The 33 entry (summed).** Per harmonic the resonant $2\pi J_n^2 D_2$ reaches the top moments
$M^3_F,M^2_T$, which carry $\zeta_n$ and diverge as $k_\parallel\!\to\!0$; the divergence cancels only after summing against $I_B$ (using $\Sigma_n J_n^2=1$). The finite `e∥e∥` density is

$$
\bigl[I(p_\perp)\bigr]_{33}
=\sum_n 2\pi J_n^2 D_2+I_B
=\sum_n 2\pi J_n^2\Bigl[\,n\Omega\,M^2_F+(\omega-n\Omega)\,p_\perp\,M^1_T\,\Bigr].
$$

All six distinct entries draw on the five parallel moments
$\{M^0_F,M^1_F,M^2_F,M^0_T,M^1_T\}$ at the slice $p_\perp$ and the six symmetric
Bessel bilinears of the regular triple $\{R_n,J_n',J_n\}$ at $z=k_\perp p_\perp/\Omega$.

#### Separable `f₀=f∥·f⊥`

Now `∂∥f₀=f∥′f⊥`, `∂⊥f₀=f∥f⊥′`, so the `p⊥`-slice pulls out of `𝒞`:

    𝒞[ p^m ∂∥f₀(·,p⊥) ](ζ_n) = f⊥(p⊥) · 𝒞[ p^m f∥′ ](ζ_n),
    𝒞[ p^m ∂⊥f₀(·,p⊥) ](ζ_n) = f⊥′(p⊥) · 𝒞[ p^m f∥ ](ζ_n).

And `X_n` **fully factors**:

    X_n = Σ_terms ( parallel 𝒞-moment of f∥ ) × ( perp P⊥-moment of f⊥ ),

**two independent 1-D primitives** — `𝒞` (§5) and `P⊥` (§6) — each computed once.
Separability would make kernel B 1-D; for A the Newberger kernel still ties `a(p∥)` to `z(p⊥)`, so separability does not factor it.

### 5.2 Relativistic Case

In $(\gamma,p_\parallel)$ coordinates

$$
I = I(\gamma)=\sum_n I_n(\gamma)+\mathbf e_\parallel\mathbf e_\parallel\,I_B(\gamma),
\qquad
I_n(\gamma)=-\frac{2\pi}{k_\parallel}\!\!\int_{|p_\parallel|<\sqrt{\gamma^2-1}}\!\!\!\!dp_\parallel\;
\frac{\mathcal U\,\boldsymbol{\mathcal T}_n}{p_\parallel-\zeta_n(\gamma)},
$$

with the non-resonant term addend ($p_\perp=\sqrt{\gamma^2-1-p_\parallel^2}$, no pole)

$$
I_B(\gamma)=2\pi\!\!\int_{|p_\parallel|<\sqrt{\gamma^2-1}}\!\!\!\!dp_\parallel\;
\Bigl(p_\parallel\,\partial_\parallel f_0-\frac{p_\parallel^2}{p_\perp}\,\partial_\perp f_0\Bigr).
$$

The inner $p_\parallel$ integral of $I_n$ is the §5 parallel Cauchy transform at the single pole
$\zeta_n(\gamma)$ (Plemelj plus growing-sheet Landau residue), over $|p_\parallel|<\sqrt{\gamma^2-1}$.

**Assembly.** Unlike the non-relativistic case there is **no factored assembly**: because
$z=(k_\perp/\Omega_0)\sqrt{\gamma^2-1-p_\parallel^2}$ couples $\gamma$ and $p_\parallel$ inside
$\boldsymbol{\mathcal T}_n$, the $p_\parallel$ integral cannot become closed moments. One simply forms the
full §3B integrand $2\pi\,\mathcal U\,\boldsymbol{\mathcal T}_n$ **pointwise** at each $(\gamma,p_\parallel)$
node and quadratures in $p_\parallel$. Two consequences: (i) the §5 summed cancellation that cures the
$m_{33}$ $k_\parallel\!\to\!0$ divergence is moot — that divergence lives in the parallel *moments*, which
are never formed here, so nothing needs cancelling. $I_B$ is the same explicit separate addend
as in the non-relativistic case;

### 5.2.1 Pushing the relativistic density with the covariant `𝒰`

The same machinery *almost* factors $I(\gamma)$. Write $\mathcal U$ in its covariant form
($\partial_\gamma$ at fixed $p_\parallel$, $\partial_{p_\parallel}$ at fixed $\gamma$)

$$
\mathcal U=\omega\,\partial_\gamma f_0+k_\parallel\,\partial_{p_\parallel}f_0,
$$

With the single-pole Cauchy transform at $\zeta_n(\gamma)$,
$\mathcal C[g]\equiv\int_{|p_\parallel|<\sqrt{\gamma^2-1}}\!g/(p_\parallel-\zeta_n(\gamma))\,dp_\parallel$,
the inner integral is:

$$
I_n(\gamma)=-\frac{2\pi}{k_\parallel}\Bigl[\,
\omega\,\mathcal C[\partial_\gamma f_0\,\boldsymbol{\mathcal T}_n]
+k_\parallel\,\mathcal C[\partial_{p_\parallel}f_0\,\boldsymbol{\mathcal T}_n]\,\Bigr].
$$

**Where it halts.** Non-relativistically the transcendental part of ${\mathcal T}_n$ (Bessel in
$z$) sat in the **outer** coordinate $p_\perp$, leaving the inner integrand polynomial in $p_\parallel$;
$\mathcal C$ then collapsed to the finite moment set $\{M^m_F,M^m_T\}$ with the Bessel weights pulled out
front. Relativistically $z=(k_\perp/\Omega_0)\sqrt{\gamma^2-1-p_\parallel^2}$ and $p_\perp$ depend on the
**inner** $p_\parallel$, so $\boldsymbol{\mathcal T}_n$ stays transcendental in $p_\parallel$. Then
The Bessel weights cannot leave the integral, and no closed $\{M^m,P_j\}$ pair survives.
The covariant $\mathcal U$ buys the clean single pole and the right derivative pair — but **not** the factorization.

Switching to $(p_\perp,p_\parallel)$ moves $z=k_\perp p_\perp/\Omega_0$ to the outer coordinate (the perp side would factor)
but then the resonance $\omega\gamma-k_\parallel p_\parallel=n\Omega_0$ becomes
quadratic in $p_\parallel$ — the relativistic resonance curve, $(\omega^2-k_\parallel^2)p_\parallel^2
-2k_\parallel n\Omega_0 p_\parallel+\omega^2(1+p_\perp^2)-n^2\Omega_0^2=0$, up to two roots per $n$ — so forfeits its single clean pole. 

Neither coordinate gives a product of independent 1-D primitives: the coupled relativistic density is irreducibly 2-D, and only $\gamma\!\to\!1$ unties $z$ from $p_\parallel$.

### 5.2.2 Edge-removing quadrature maps

Both $I_n(\gamma)$ and $I_B(\gamma)$ integrate the momentum disk $|p_\parallel|<u_{\max}$, $u_{\max}\equiv\sqrt{\gamma^2-1}$,
with $p_\perp=\sqrt{u_{\max}^2-p_\parallel^2}$, and the outer $\gamma\in[1,\gamma_{\max}]$. Two square-root edges
make the integrand non-smooth on the boundary and stall any quadrature that resolves them directly:

- **Rim** $p_\perp\to0$ (at $|p_\parallel|=u_{\max}$): the numerator carries an explicit $1/p_\perp$ —
  $\mathcal U\supset(\omega\gamma-k_\parallel p_\parallel)\,p_\perp^{-1}\partial_\perp f_0$ and
  $I_B\supset p_\parallel^2 p_\perp^{-1}\partial_\perp f_0$ — and $p_\perp(p_\parallel)=\sqrt{u_{\max}^2-p_\parallel^2}$
  has a vertical tangent there.
- **Floor** $\gamma\to1$: $u_{\max}=\sqrt{\gamma^2-1}\sim\sqrt{2(\gamma-1)}$ collapses the disk with a
  $\sqrt{\gamma-1}$ edge.

Two substitutions remove both, mapping the disk to a fixed square $(q,\theta)\in[0,1]\times[-\tfrac\pi2,\tfrac\pi2]$:

$$
p_\parallel=u_{\max}\sin\theta,\quad p_\perp=u_{\max}\cos\theta,\quad dp_\parallel=u_{\max}\cos\theta\,d\theta\;(=p_\perp\,d\theta);
\qquad
\gamma=1+(\gamma_{\max}-1)\,q^2,\quad d\gamma=2(\gamma_{\max}-1)\,q\,dq.
$$

The inner Jacobian $u_{\max}\cos\theta=p_\perp$ **cancels the $1/p_\perp$** exactly and renders $p_\perp$ the
analytic $u_{\max}\cos\theta$; the outer $q^2$ map makes $u_{\max}\propto q$ and $d\gamma\propto q\,dq$, so the
integrand $\times$ Jacobian vanishes smoothly as $q\to0$. A relativistic contribution becomes

$$
\int_1^{\gamma_{\max}}\!\!\!d\gamma\!\!\int_{-u_{\max}}^{u_{\max}}\!\!\!\!dp_\parallel\,(\cdots)
=\int_0^1\!\!\!dq\!\!\int_{-\pi/2}^{\pi/2}\!\!\!\!d\theta\;
\underbrace{2(\gamma_{\max}-1)\,q\;u_{\max}\cos\theta}_{\text{Jacobian}}\;(\cdots),
$$

with the integrand now smooth on the closed box. The analytic Landau term
$\propto\log\frac{u_{\max}-\zeta_n}{-u_{\max}-\zeta_n}$ (the pole part of $\mathcal C$, §5.3) is left untouched —
only the regularized remainder is quadratured. Fixed Gauss–Legendre on the box then converges geometrically;
on the raw disk the boundary edges cap it near $\sim\!10^{-4}$.

### 5.3 Parallel Cauchy transform and resonance moments

The parallel primitive is the Landau-causal Cauchy transform over `supp g` (finite per cell, infinite for a Gaussian; `Im ζ>0` physical):

$$
\mathcal C[g](\zeta) \equiv \int \frac{g(p)}{p-\zeta}\,dp,
\qquad \zeta_n=\frac{\omega-n\Omega}{k_\parallel}.
$$

The surrounding derivation uses this same `𝒞`. “Cauchy transform” is the precise
off-real-axis name; its real-axis boundary value is a Hilbert transform plus the Plemelj
residue, so the code field `hilbert` is a (retained) misnomer for the off-axis object. The actual resonance moments are the single
indexed functional

$$
\mathcal C_m[g](\zeta_n)
\equiv -\frac1{k_\parallel}\mathcal C[p^m g](\zeta_n)
=\int\frac{p^m g(p)}{\omega-k_\parallel p-n\Omega}\,dp.
$$

For a separable distribution, the two inputs are `g=f∥` and `g=f∥′`:

$$
M_F^m\equiv\mathcal C_m[f_\parallel],\qquad
M_T^m\equiv\mathcal C_m[f_\parallel'].
$$

Thus `F` and `T` are implementation-facing labels, not distinct transforms. No universal
symbol exists for a distribution-agnostic Cauchy moment; $\mathcal C_m[g]$ states its
meaning directly. Reserve $Z_m$ for Maxwellian-only moments. Use `M_F^m,M_T^m` only where
matching the assembler or code fields helps. The `1/k∥` and sign come from
`ω−k∥p−nΩ = −k∥(p−ζ)`. Both are evaluations of the *same* $\mathcal C$, dispatched on
how `g` is represented:

| `g` | $\mathcal C[g]$ |
|---|---|
| normalized Gaussian `e^{−p²}/√π` | `Z(ζ)`, plasma dispersion function |
| unnormalized Gaussian `e^{−p²}` | `√π Z(ζ)` |
| piecewise polynomial | per cell `∫q dp + P(ζ)·log((p_{i+1}−ζ)/(p_i−ζ))` |
| arbitrary analytic | Plemelj split (below) |

**Piecewise-poly cell.** Synthetic-divide `P(p)=q(p)(p−ζ)+P(ζ)`; then
`∫P/(p−ζ)dp = ∫q dp + P(ζ)·log((p_{i+1}−ζ)/(p_i−ζ))`. **Branch-cut invariant**: one
complex `log` of the *ratio* (not a difference of logs) keeps the continuation
single-valued as `Im ζ→0`. Landau continuation to the growing sheet (`Im ζ<0`, `Re ζ` in
cell) adds `2πi·p(ζ)`.

**Arbitrary analytic `g`.** Plemelj split with the removable singularity pulled out:

    𝒞[g](ζ) = ∫ (g(p)−g(ζ))/(p−ζ) dp + g(ζ)·log((b−ζ)/(a−ζ)) [ + 2πi·g(ζ) ]

The first integrand is regular at `p=ζ` ⇒ plain adaptive quadrature; the `log` ratio
carries the branch cut; the `Im ζ→0⁺` limit supplies the Plemelj `+iπ g(ζ)`, and the
explicit `2πi g(ζ)` is the lower-half Landau term.

### 5.4 Perpendicular primitive `P⊥`

The perpendicular primitive collects the `p⊥` integrals (`2π` = azimuthal factor of
`d³p`). The Bessel argument `z=k⊥p⊥/Ω₀` runs with `p⊥` (`Jn′≡∂Jn`). Index the three
Bessel bilinears by `j=0,1,2`, `Wⱼ∈{Jₙ², JₙJₙ′, Jₙ′²}`; the `p⊥`-power is fixed by `j` and
the density (`T_n` entries, §3B). The **perp moments**:

$$
P_j \equiv 2\pi\!\int W_j(z)\,f_\perp\,p_\perp^{\,j+1}\,dp_\perp,\qquad
P_j^\partial \equiv 2\pi\!\int W_j(z)\,f_\perp'\,p_\perp^{\,j}\,dp_\perp,
$$

(`∂` superscript = the `f⊥′` slice). The six, with code fields:

| `j` | `Wⱼ` | `Pⱼ` (code) | `Pⱼ^∂` (code) |
|---|---|---|---|
| 0 | `Jₙ²`   | `JF`  = 2π∫Jₙ²f⊥ p⊥ dp⊥     | `J∂F`   = 2π∫Jₙ²f⊥′ dp⊥     |
| 1 | `JₙJₙ′` | `JdJF`= 2π∫JₙJₙ′f⊥ p⊥² dp⊥  | `JdJ∂F` = 2π∫JₙJₙ′f⊥′ p⊥ dp⊥ |
| 2 | `Jₙ′²`  | `∂J²F`= 2π∫Jₙ′²f⊥ p⊥³ dp⊥   | `∂J²∂F` = 2π∫Jₙ′²f⊥′ p⊥² dp⊥ |

Dispatched on `f⊥`:

| `f⊥` | moments |
|---|---|
| Gaussian | `Γ_n(λ)=Iₙ(λ)e^{−λ}`, `λ=(k⊥ p_th⊥/Ω₀)²/2`, plus recurrences |
| piecewise polynomial | per-cell Bessel-product power series (Schläfli ₂F₃) `∫pᵈ⁺¹ J J` |
| arbitrary analytic | direct adaptive quadrature of the six integrals |

The Gaussian closed forms `{P₀=Γ_n, P₀^∂=−2Γ_n/p_th², P₁^∂=−(k⊥/Ω₀)Γ_n′, …}` are the
normalized `f⊥` evaluations of exactly these integrals — closed per harmonic, but the
**sum over `n`** of `Γ_n` is still truncated (Gaussian is evaluator B).


## 7. Worked example: drifting bi-Maxwellian closes both primitives

Take the normalized drifting bi-Maxwellian `f₀=f∥·f⊥`,

$$
f_\|(p)=\frac{e^{-(p-v_d)^2/p_{\mathrm{th}\|}^2}}{\sqrt\pi\,p_{\mathrm{th}\|}},\qquad
f_\perp(p)=\frac{e^{-p^2/p_{\mathrm{th}\perp}^2}}{\pi\,p_{\mathrm{th}\perp}^2}.
$$

Two master integrals close the primitives, and in each the derivative form needs *no* new
integral (since `f'=−2(p−v_d)f/p_th²`).

**Parallel `𝒞` → plasma dispersion function `Z`.** Rescale `u=(p−v_d)/p_th∥`; the pole
maps to

$$
\xi_n=\frac{\zeta_n-v_d}{p_{\mathrm{th}\|}}=\frac{\omega-k_\| v_d-n\Omega_0}{k_\| p_{\mathrm{th}\|}},
\qquad
\mathcal C[f_\|](\zeta_n)=\frac1{p_{\mathrm{th}\|}}\,Z(\xi_n).
$$

The derivative collapses onto the same family via `f∥′=−2(p−v_d)f∥/p_th∥²`:

$$
\mathcal C[f_\|'](\zeta_n)=-\frac2{p_{\mathrm{th}\|}^2}\bigl[1+\xi_n Z(\xi_n)\bigr]=\frac1{p_{\mathrm{th}\|}^2}Z'(\xi_n),
$$

giving

$$
M^0_F=-\frac{Z(\xi_n)}{k_\| p_{\mathrm{th}\|}},\qquad
M^0_T=-\frac{Z'(\xi_n)}{k_\| p_{\mathrm{th}\|}^2}.
$$

**Perpendicular `P⊥` → ring sum `Γ_n`.** Weber's second exponential integral,

$$
\int_0^\infty e^{-p^2/p_{\mathrm{th}\perp}^2}J_n(\beta p)^2\,p\,dp
=\frac{p_{\mathrm{th}\perp}^2}{2}\,e^{-\lambda}I_n(\lambda),
\qquad \beta=\frac{k_\perp}{\Omega_0},\quad
\lambda\equiv\frac{(\beta\,p_{\mathrm{th}\perp})^2}{2},
$$

collapses the base moment to the modified-Bessel ring sum (`z=k⊥p/Ω₀`):

$$
P_0=2\pi\!\int_0^\infty J_n(z)^2 f_\perp\,p\,dp=e^{-\lambda}I_n(\lambda)\equiv\Gamma_n(\lambda).
$$

Again the derivative needs no new integral — `f⊥′=−2p f⊥/p_th⊥²` gives

$$
P_0^\partial=2\pi\!\int_0^\infty J_n(z)^2 f_\perp'\,dp=-\frac{2}{p_{\mathrm{th}\perp}^2}\,P_0=-\frac{2\,\Gamma_n(\lambda)}{p_{\mathrm{th}\perp}^2},
$$

and the mixed combo from differentiating Weber w.r.t. `β`: `P₁^∂=−(k⊥/Ω₀)Γ_n′`, with
`Γ_n′(λ)=½(Γ_{n−1}+Γ_{n+1})−Γ_n`.