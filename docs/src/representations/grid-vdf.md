# GridVDF

`GridVDF` projects tabulated gyrotropic data
`F[i,j] = f₀(v⊥ᵢ,v∥ⱼ)` onto a compactly supported callable distribution.

## Projection

`fit_grid(method, vperp, vpara, F)` returns a projection implementing:

- `projection(v⊥, v∥)`
- `domain(projection)`
- optionally specialized `_grad2` and `moment`

`NonnegBSpline` is the default: an adaptive, nonnegative tensor B-spline fit.
`BicubicHermite` provides local C¹ interpolation without a positivity guarantee.
Both produce `TensorSplineFit`, stored as local power coefficients on each cell.

`moment` supplies density and mean squared perpendicular momentum. Tensor fits integrate
these moments exactly; other projections use quadrature. Susceptibility is divided
by cached density, so input amplitude is irrelevant.

## Tensor fast path

Non-relativistic `HarmonicSum` evaluation exploits the cell polynomials. For
harmonic `n`, the parallel pole is

```math
\zeta_n = \frac{\omega-n\Omega}{k_\parallel}.
```

For a parallel-cell polynomial `P`,

```math
\int_l^r \frac{P(u)}{u-\zeta}\,du
= \int_l^r Q(u)\,du
+ P(\zeta)\log\frac{r-\zeta}{l-\zeta},
```

where `P(u) = Q(u)(u-\zeta) + P(ζ)`. `_cell_hilbert` evaluates this form and
adds the oriented Landau jump when the pole crosses the cell.

Within one perpendicular cell, every parallel slice is polynomial in the local
coordinate `t`. Linearity of the Hilbert integral therefore makes all five
parallel moments polynomials in `t`. Their coefficients and shared cell logarithms
are computed once per harmonic and perpendicular cell, not once per quadrature
node.

The remaining perpendicular integral evaluates those polynomials with Bessel
weights using adaptive Gauss–Kronrod quadrature. Cells are split at the Bessel
oscillation scale. At `k∥ = 0`, plain cell moments replace the Hilbert integral.

Relativistic evaluation and non-tensor projections use the generic `CoupledVDF` path.
