
Relativistic f₀ works better in $(\gamma,p_\parallel)$ coordinates.

The denominator is nonlinear in `p∥` through `γ`, but at fixed `γ` it linearizes:

    ω − k∥v∥ − nΩ = (ωγ − k∥p∥ − nΩ₀)/γ = −(k∥/γ)(p∥ − ζ_n(γ)),   ζ_n(γ) = (ωγ − nΩ₀)/k∥,

with a clean rational pole in `p∥`:

    χ_n = (Π²/ω²)·(−2π/k∥) ∫ dγ ∫_{|p∥|<√(γ²−1)} dp∥ · 𝒰 𝓣_n /(p∥ − ζ_n(γ)).

This straightens A's `sin πa=0` resonance curve into the line `p∥=ζ_n(γ)`. What stays
coupled is only `z=(k⊥/Ω₀)√(γ²−1−p∥²)`, so the `p∥` integral is the analytic `𝒞` branch over a finite interval.


The outer-coordinate density `I` is given by

$$
I = I(\gamma)=\sum_n I_n(\gamma)+\mathbf e_\parallel\mathbf e_\parallel\,I_B(\gamma),
\\
I_n(\gamma)=-\frac{2\pi}{k_\parallel} \int dp_\parallel
\frac{\mathcal U\, {\mathcal T}_n}{p_\parallel-\zeta_n(\gamma)}, \quad
I_B(\gamma)=2\pi\int dp_\parallel
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
