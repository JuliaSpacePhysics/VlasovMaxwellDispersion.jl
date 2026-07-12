# Relativistic path: (p⊥,p∥) slicing

The relativistic resonant denominator, in momentum variables $p/mc$ with $\gamma=\sqrt{1+p_\perp^2+p_\parallel^2}$,

$$
\omega - k_\parallel v_\parallel - \frac{n\Omega_0}{\gamma} = \frac{D_n}{\gamma},\qquad
D_n = \omega\gamma - k_\parallel p_\parallel - n\Omega_0,
$$

gives per harmonic the term the **`HarmonicSum`** backend sums, $\chi=\sum_n\chi_n$,

$$
\chi_n \propto \int_0^{p_{\perp\max}}\!\!dp_\perp \int_{-P}^{P}\!\!dp_\parallel\;
\frac{\mathcal U\,\boldsymbol{\mathcal T}_n\;p_\perp/\gamma}{D_n},
$$

with the covariant numerator $\mathcal U = k_\parallel\partial_\parallel f_0 + (\omega\gamma-k_\parallel p_\parallel)\,p_\perp^{-1}\partial_\perp f_0$
and $\boldsymbol{\mathcal T}_n(z,p_\parallel,p_\perp)$ the harmonic Bessel tensor.

The Bessel argument $z = k_\perp p_\perp/\Omega_0$ is **independent of the inner
$p_\parallel$**, so $J$'s are computed once per slice. 
The non-resonant $\mathbf e_\parallel\mathbf e_\parallel$ term adds a smooth integral:
$I_B = 2\pi\iint (p_\perp p_\parallel \partial_\parallel f_0 - p_\parallel^2\,\partial_\perp f_0)/\gamma$.

The Qin's formulation uses the complex-order $\sigma$-quartet
$\boldsymbol{\mathcal T}(a,z)$, $a=(\omega\gamma-k_\parallel p_\parallel)/\Omega_0$.

## Pole structure: rationalization

$D_n$ is not rational in $p_\parallel$, but

$$
D_n\tilde D_n = \omega^2\gamma^2 - (k_\parallel p_\parallel + n\Omega_0)^2
= A\,p_\parallel^2 + B\,p_\parallel + C,\qquad \tilde D_n = \omega\gamma + k_\parallel p_\parallel + n\Omega_0,
$$

with $A=\omega^2-k_\parallel^2$, $B=-2k_\parallel n\Omega_0$, $C=\omega^2 m_\perp^2 - n^2\Omega_0^2$,
$m_\perp^2 \equiv 1+p_\perp^2$. So $1/D_n = \tilde D_n/[A(p_\parallel-p_+)(p_\parallel-p_-)]$: two
explicit simple poles, and the **squaring ghost** — the quadratic root that solves
$\tilde D_n=0$ — carries an identically-zero residue.
No root classification is ever needed. 
Residues are $r_\pm = \mp g(p_\pm)\tilde D_n(p_\pm)/\sqrt{B^2-4AC}$; near-axis poles are peeled by the standard Plemelj subtraction, far ones left to plain quadrature.

## The Landau rule: poles cross the axis only at Im ω = 0

For $\omega=\omega_r+i\nu$, $B$ is real and $\operatorname{Im}A = 2\omega_r\nu$,
$\operatorname{Im}C = 2\omega_r\nu\,m_\perp^2$; a real root would need
$\operatorname{Im}(Ap^2+Bp+C) = 2\omega_r\nu\,(p^2+m_\perp^2)=0$, i.e. $p^2=-m_\perp^2$ —
impossible. **Poles touch the real $p_\parallel$ axis only at $\nu=0$**.
Hence the complete continuation bookkeeping (continuation defined from $\nu\to+\infty$,
$k_\parallel>0$ convention):

- $\nu>0$: straight integral, no residue terms;
- $\nu=0$: boundary value from each pole's home side, given exactly by the local slope $dp/d\omega = \gamma^2/(k_\parallel\gamma-\omega p)$;
- $\nu<0$: $+2\pi i\,r$ for in-range poles found below the axis. This is the exact
  continuation **iff subluminal** ($\omega_r^2<k_\parallel^2$): there
  $\operatorname{disc} = B^2-4AC = 4\omega^2[\,n^2\Omega_0^2+(k_\parallel^2-\omega^2)m_\perp^2\,]$
  is uniformly nonzero, the true pole has home side $+$ for all $(n,p_\perp)$
  ($k_\parallel\gamma-\omega p>0$ on the physical branch) and crosses downward, while the
  ghost crosses upward with null coefficient.

Support endpoints $|p_\parallel|=P$ sit where $f_0\approx 0$, so pole–endpoint collisions
carry negligible coefficients: no endpoint corrections.

## Why damped superluminal is special

Superluminally ($\omega_r^2>k_\parallel^2$) the resonance curve is an ellipse and
$\operatorname{disc}$ vanishes at its apex,

$$
m_\perp^{*2}(\omega) = \frac{n^2\Omega_0^2}{\omega^2-k_\parallel^2},
$$

real and inside the support for harmonics in the resonant band
$k_\parallel^2<\omega^2<k_\parallel^2+n^2\Omega_0^2$. There the two poles coalesce — a
**pinch** of the $p_\parallel$ contour — and the slice function $H_n(p_\perp,\omega)$ has a
branch point at $p_\perp^*(\omega)$. Its migration under $\omega\to\omega+i\nu$ follows
$\operatorname{Im}m_\perp^{*2} = -2\omega_r\nu\,n^2\Omega_0^2/|\omega^2-k_\parallel^2|^2$:

- $\nu>0$: the branch point sits **below** the real $p_\perp$ path — straight integral is
  analytic, equals $\chi_n$. At $\nu=0$
  the apex is an integrable $1/\sqrt{\,}$ kink on the path;
- $\nu<0$: it **crosses to above the path**, so the straight real-sliced integral is no
  longer the analytic continuation — the missing piece is the cut discontinuity around
  $p_\perp^*(\omega)$.

This is a geometric invariant, not a coordinate artifact: the apex lies at
$\gamma = n\Omega_0\omega/(\omega^2-k_\parallel^2)$ — the vertex of the same quadratic
$\gamma^2(\omega^2-k_\parallel^2) - 2\gamma n\Omega_0\omega + n^2\Omega_0^2+k_\parallel^2=0$
whose roots were the rim crossings of the retired $(\gamma,p_\parallel)$ slicing. Any real
slicing of the 2-D momentum integral meets this tangency for damped superluminal $\omega$
and needs either a cut correction or a complex path.

Rather than re-grow cut-correction machinery for this physically thin regime (heavily
damped EM branches; the Swanson closed form is not the continuation there either and
warns), damped superluminal $\omega$ is **not supported**: `contribution` warns and
returns the straight integral, which is off by the O(1) missed crossings. To reach such
roots, evaluate on $\operatorname{Im}\omega\ge0$ — where the path is exact — and continue
externally (least-squares polynomial in $\operatorname{Im}\omega$; recipe and
cross-grid certification in `test-relativistic.jl`). Building this into the evaluator
was tried and dropped: correct, but 12 samples per call made root polishing ~200× slower
and the fit noise capped the achievable residual.

Subluminally the apex roots give $p_\perp^{*2}\le -1$: imaginary, at distance $\ge1$ from
the path — uniformly harmless, which is why the subluminal rule above is exact at any
damping depth reached before the poles approach the $\gamma$ branch points at
$p_\parallel=\pm i\,m_\perp$ ($|\operatorname{Im}p_\parallel|\sim|\nu|\gamma/|k_\parallel-\omega v_\parallel|\lesssim 1$;
validated to $\operatorname{Im}\tilde\omega=-0.15$ at $\mu=2$, where the old rim-corrected
path had drifted to $2.6\times10^{-4}$).

Caveat: for **tabulated** f₀ (GridVDF) the damped-side residues probe the spline's
analytic continuation at $|\operatorname{Im}p_\parallel|\sim|\nu|\gamma/|k_\parallel-\omega v_\parallel|$
— several cell-widths off-axis by $\nu\sim-0.02$, where piecewise-polynomial
continuation amplifies fit noise (observed floor $\sim3\times10^{-2}$ at $\mu=2$ on a
61×121 grid vs $10^{-5}$ for analytic f₀). Finer grids push it down.
