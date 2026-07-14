# Compressing general coupled VDF with Adaptive Cross Approximation (ACA)

This approach approximates the VDF or its gradients on a real (q,u) grid by low-rank crossed slices (ACA/CUR matrix approximation), reusable over many k. Landau transforms are applied analytically to separated parallel factors.

It trades the exact path's 2-D adaptive quadrature per $\omega$ for a few Landau integrals.

## Notation

| Symbol | Meaning |
| --- | --- |
| $\Omega,\ \Pi^2$ | Signed gyrofrequency $\Omega_s/\Omega_\text{ref}$ and $(\omega_{ps}/\Omega_\text{ref})^2$ |
| $(v,u)$ | Non-relativistic perpendicular / parallel velocity |
| $f_0(v,u)$ | Gyrotropic equilibrium; $N=\int 2\pi v\,f_0\,dv\,du$ |
| $a,\ z$ | Bessel scale $a=k_\perp/\Omega$, argument $z=av$ |
| $n$ | Cyclotron harmonic, summed over $-n_{\max}..n_{\max}$ |
| $D_n(u)$ | Resonant denominator $\omega-n\Omega-k_\parallel u = -k_\parallel(u-\zeta_n)$ |
| $\zeta_n$ | Parallel pole $(\omega-n\Omega)/k_\parallel$ |
| $\mathcal L,\ \sigma$ | Landau contour; causal side $\sigma\,\mathrm{Im}\,\zeta_n>0$, $\sigma=\operatorname{sign}k_\parallel$ |

## 1. Starting point

The susceptibility is a cyclotron-harmonic sum whose integrand is bilinear: 
a $u$-independent Bessel tensor $\mathbf K_n(v)$ times a gradient
numerator $\mathcal G_n(v,u)$ built from $\partial_v f_0,\ \partial_u f_0$,

```math
\chi_s=\frac{\Pi^2}{N\omega^2}\sum_n\int 2\pi v\,dv\int
\frac{\mathbf K_n(v)\,\mathcal G_n(v,u)}{D_n(u)}\,du .
```

Only $\mathbf K_n(v)$ survives the $u$-integral, so it can be precomputed; $\mathcal G_n(v,u)$ is what
the parallel integral acts on each $\omega$. ($\mathbf K_n$ is defined in §4.)

## 2. Rank-$R$ separation

The cross writes $f_0$ as $R$ separable terms whose parallel factors $b_s(u)=f_0(v_s,u)$ are **true
slices** of $f_0$ — analytic, hence continuable off the real axis:

```math
f_0(v,u)\approx\sum_{s=1}^R \tilde a_s(v)\,b_s(u),\qquad
\partial_v f_0=\sum_s \tilde a_s'(v)\,b_s(u),\qquad
\partial_u f_0=\sum_s \tilde a_s(v)\,b_s'(u).
```

Now $v$ and $u$ separate, so every $(n,s)$ term splits into a perpendicular tensor integral
(precomputable, §4) and a parallel moment integral (per $\omega$, §3).

## 3. Parallel moments

The $u$-integral of $\mathcal G_n$ collapses to two moment families — one from the slice $b_s$, one
from its derivative $b_s'$:

```math
\Phi_m=\int_{\mathcal L}\frac{u^m b_s(u)}{D_n(u)}\,du\ \ (m=0,1,2),\qquad
\Psi_m=\int_{\mathcal L}\frac{u^m b_s'(u)}{D_n(u)}\,du\ \ (m=0,1).
```

Only $m=0$ is a genuine transform. Using $u^m=\zeta_n u^{m-1}+u^{m-1}(u-\zeta_n)$ and the
$\omega$-independent raw moments $I_m=\int u^m b_s\,du$, $J_m=\int u^m b_s'\,du$,

```math
\Phi_m=\zeta_n\Phi_{m-1}-\frac{I_{m-1}}{k_\parallel},\qquad
\Psi_m=\zeta_n\Psi_{m-1}-\frac{J_{m-1}}{k_\parallel}.
```

When $k_\parallel=0$ the denominator $D_n=\omega-n\Omega$ is constant, so no transform is needed:
$\Phi_m=I_m/(\omega-n\Omega)$, $\Psi_m=J_m/(\omega-n\Omega)$.

### The base transform

$\Phi_0=-\mathcal C[b_s]/k_\parallel$ and $\Psi_0=-\mathcal C[b_s']/k_\parallel$ are the Landau-Cauchy
transform of a slice $g\in\{b_s,b_s'\}$ along the contour $\mathcal L$ over the support $[l,h]$ —
the principal-value integral plus the residue of any pole that damping has dragged onto the causal
side:

```math
\mathcal C[g](\zeta_n)=\int_{\mathcal L}\frac{g(u)}{u-\zeta_n}\,du
=\;\mathrm{p.v.}\!\int_l^h\frac{g(u)}{u-\zeta_n}\,du
\;+\;\sigma\,2\pi i\,g(\zeta_n)\,\Big|_{\text{crossed}},
```

with *crossed* meaning $\mathrm{Re}\,\zeta_n\in(l,h)$ and $\sigma\,\mathrm{Im}\,\zeta_n<0$. Because $g$
is a true slice of $f_0$, $g(\zeta_n)$ — and thus the Landau residue — is exact; a fitted (spline /
SVD) surrogate has no value off the real axis and returns wrong growth rates. This is the whole
reason the surrogate reaches damped modes.

## 4. Perpendicular tensors (precomputed once per $\mathbf k$)

With $z=av$,

```math
\mathbf b_n(v)=\bigl(v R_n,\ v J_n',\ J_n\bigr),\qquad
R_n=\tfrac12(J_{n-1}+J_{n+1})=\tfrac{n}{z}J_n,\qquad
J_n'=\tfrac12(J_{n-1}-J_{n+1}),
```

and $\mathbf K_n(v)=\mathbf b_n\mathbf b_n^{\top}$ (6 distinct entries). Pairing $\mathbf K_n$ with the
$\partial_v$ slice $\tilde a_s'$ and the $b_s'$-partnered factor $v\tilde a_s$ gives two
$\omega$-independent tensors:

```math
\mathbf P^{\partial}_{n,s}=2\pi\int \tilde a_s'(v)\,\mathbf K_n(v)\,dv,\qquad
\mathbf P^{F}_{n,s}=2\pi\int v\,\tilde a_s(v)\,\mathbf K_n(v)\,dv .
```

## 5. Assemble the tensor

Parallel weights combine the two slices at each order:

```math
w_F^0=\omega\Phi_0-k_\parallel\Phi_1,\quad
w_F^1=\omega\Phi_1-k_\parallel\Phi_2,\quad
w_T^0=k_\parallel\Psi_0,\quad
w_T^1=k_\parallel\Psi_1 .
```

Each independent entry is a perp-tensor component times its weight, summed over $(n,s)$:

```math
\begin{aligned}
\chi_{xx}&=P^{\partial}_{11}w_F^0+P^{F}_{11}w_T^0, &
\chi_{xy}&=i\bigl(P^{\partial}_{12}w_F^0+P^{F}_{12}w_T^0\bigr), &
\chi_{yy}&=P^{\partial}_{22}w_F^0+P^{F}_{22}w_T^0,\\
\chi_{xz}&=P^{\partial}_{13}w_F^1+P^{F}_{13}w_T^1, &
\chi_{yz}&=i\bigl(P^{\partial}_{23}w_F^1+P^{F}_{23}w_T^1\bigr), &
\chi_{zz}&=n\Omega\,P^{\partial}_{33}\Phi_2+(\omega-n\Omega)P^{F}_{33}\Psi_1 .
\end{aligned}
```

With the prefactor, the gyrotropic tensor is

```math
\chi_s=\frac{\Pi^2}{N\omega^2}
\begin{pmatrix}
\chi_{xx} & \chi_{xy} & \chi_{xz}\\
-\chi_{xy} & \chi_{yy} & -\chi_{yz}\\
\chi_{xz} & \chi_{yz} & \chi_{zz}
\end{pmatrix},
```

the off-diagonal sign pattern being the gyrotropic symmetry of $\chi_s$.

## Implementation

The mathematics above is exact; the following are numerical choices in the code.

**Evaluating Landau integral $\mathcal C$ — far/near split.** $\mathcal C$ is computed by one of two schemes
selected by the pole distance $|\zeta_n|$ against $\theta U$ ($U=\max(|l|,|h|)$, $\theta=2$). Both
return the same transform (each adds the crossing residue of §3); the split is purely for cost and
conditioning, and is where the speedup lives.

- **Far, $|\zeta_n|>\theta U$:** $u/\zeta_n$ is small on $[l,h]$, so
  $\tfrac1{u-\zeta_n}=-\sum_{p\ge0}u^p/\zeta_n^{p+1}$ gives a Neumann series in the
  $\omega$-independent moments $\int_l^h u^p g\,du$ (built once) —
  $\ \mathrm{p.v.}\!\int=-\sum_{p\ge0}\zeta_n^{-(p+1)}\int_l^h u^p g\,du$, truncated at $p=40$. No logs,
  no cancellation; a far harmonic costs a dot product. At quasi-perpendicular propagation nearly every
  harmonic is far — hence the win.
- **Near, $|\zeta_n|\le\theta U$:** pole subtraction on fixed Gauss nodes,
  $\ \mathrm{p.v.}\!\int=\int_l^h\frac{g(u)-g(\zeta_n)}{u-\zeta_n}\,du+g(\zeta_n)\ln\frac{h-\zeta_n}{l-\zeta_n}$.
  If $g(\zeta_n)$ grows large off-axis the subtraction cancels catastrophically, so the direct
  unpeeled sum is used instead.

**Cost.** Everything in §4 precomputes. Each $\omega$ costs, per $(n,s)$: the two base transforms
$\Phi_0,\Psi_0$, three scalar recursions, and one tensor contraction — $\mathcal O(R\,n_{\max})$ total,
against the exact [`CoupledVDF`](@ref) path's 2-D adaptive quadrature at every $\omega$.
