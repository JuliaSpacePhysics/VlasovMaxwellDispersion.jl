# GridVDF

How a tabulated gyrotropic `f₀` on a velocity grid becomes a magnetized
susceptibility `χ(ω,k)`. Non-relativistic, oblique (`k⊥≠0`), harmonic-sum closure. 
Notation: `Ω` signed gyrofrequency, `Π²` plasma-frequency ratio, `z = k⊥v⊥/Ω`.

## 1. Representation: grid → tensor spline / basis projection

Fit a tabulated `f₀[i,j] = f₀(v⊥ᵢ, v∥ⱼ)` to a tensor spline on cells `[uᵢ,uᵢ₊₁]×[wⱼ,wⱼ₊₁]`.

Default: **non-negative least-squares tensor B-spline** — adaptive knots to a
relative-error tolerance, positivity by construction, analytic derivatives. Naïve
per-cell interpolation gives negative f₀ and noisy ∂f₀ near the pole and is
rejected.

The fit is renormalized so `∫d³p f₀ = 1` (closed form over the cells, `_fit_d3p`).

## 2. Susceptibility as parallel × perpendicular moments


Each cyclotron tensor `χ_n` is a **bilinear form** in five parallel moments and six
Bessel-weighted perpendicular moments:

**Parallel moments** (the only place the Landau pole lives) — Hilbert transforms of the `v⊥`-slice of `∇f₀`:

    z^p_F(v⊥) = −(1/k∥) ∫ v∥^p ∂⊥f₀(v∥,v⊥) /(v∥ − ζ_n) dv∥,   p = 0,1,2
    z^p_T(v⊥) = −(1/k∥) ∫ v∥^p ∂∥f₀(v∥,v⊥) /(v∥ − ζ_n) dv∥,   p = 0,1

**Perpendicular weights** carry the Bessel content, `z = k⊥v⊥/Ω`:

    2π v⊥ Jₙ², 2π Jₙ², 2π v⊥ JₙJₙ′, 2π v⊥² JₙJₙ′, 2π v⊥² Jₙ′², 2π v⊥³ Jₙ′².

    χ_n = ∫₀^∞ dv⊥  𝔅( {z^p_F(v⊥), z^p_T(v⊥)},  {Bessel weights}(v⊥) ),

with `𝔅` the 3×3 entrywise sum-of-products. `converge` truncates the `n`-sum at `nmax ≈ k⊥ρ` (`nmax_bessel`).

## 3. Exact parallel primitive (piecewise polynomial)

Because `f₀` is piecewise-polynomial, each `z^p` closes **exactly** per parallel cell. For a cell polynomial `P(v∥)`, synthetic-divide by `(v∥−ζ)`,
`P = q·(v∥−ζ) + P(ζ)`:

    ∫_{uᵢ}^{uᵢ₊₁} P(v∥)/(v∥−ζ) dv∥ = ∫q dv∥ + P(ζ)·log((uᵢ₊₁−ζ)/(uᵢ−ζ))   [+ 2πi·P(ζ)]

The single complex `log` of the ratio keeps the continuation single-valued as
`Im ζ→0`; the `+2πi·P(ζ)` is the Landau term for a lower-half pole inside the cell
(`cell_hilbert_landau`). So

    z^p_F(v⊥) = −(1/k∥) Σ_i  cellH( v∥^p · ∂⊥f₀-slice ;  uᵢ, uᵢ₊₁, ζ_n ),

and likewise `z^p_T` with `∂∥f₀`.

## 4. Optimization: the parallel moment is a polynomial in v⊥

The key structural fact that makes this fast. Fix the perp cell `j` and let
`t = v⊥ − wⱼ`. The slice coefficients are polynomials in `t`:

    ∂⊥f₀ slice, coeff of s∥^{A-1} = Σ_{B≥2}(B-1) c[i,j,A,B] t^{B-2}   (deg 2 in t)
    ∂∥f₀ slice, coeff of s∥^{A-2} = (A-1) Σ_B c[i,j,A,B] t^{B-1}      (deg 3 in t)

`cellH` is **linear** in the cell polynomial's coefficients, so each parallel moment is itself a polynomial in `t`:

    z^p_F(t) = Σ_{b=0}^{2} μ^{p,b}_F · t^b,     z^p_T(t) = Σ_{b=0}^{3} μ^{p,b}_T · t^b,
    μ^{p,b} = −(1/k∥) Σ_i cellH( v∥^p · [t^b-component of the slice] ; uᵢ, uᵢ₊₁, ζ_n ).

Two precomputations, **once per (harmonic n, perp cell j)** instead of per quadrature node:

1. the moment coefficients `μ^{p,b}` (a handful of `cellH` per parallel cell);
2. the per-parallel-cell `log((uᵢ₊₁−ζ_n)/(uᵢ−ζ_n))` and Landau flag — these depend only on `(i, ζ_n)`, not on `v⊥` or the moment, so they are shared across all moments, `t`-powers and perp cells of the harmonic (`_cellH` takes the log as an argument and reuses the Horner remainder `pζ = P(ζ)` for the Landau term).

The remaining per-harmonic work is a smooth Gauss–Kronrod over each perp cell whose integrand only **evaluates the cubics** `evalpoly(t, μ)` and the Bessel weights:

    χ_n = Σ_j ∫_{wⱼ}^{wⱼ₊₁} dv⊥  𝔅( {evalpoly(t,μ_F), evalpoly(t,μ_T)}, {Bessel}(v⊥) ).

**Cost.** Per harmonic the `cellH` count drops from `O(N_⊥nodes · N_∥cells)` (the slice
re-summed at every `v⊥` node) to `O(N_∥cells)` precompute; the `N_⊥nodes` inner
evaluations become cubic + Bessel only. Measured ~**3–5× per `contribution`**
(`exact = exact`, agrees with the independent `CoupledVDF` path to ~1e-8), which is the dominant per-evaluation cost of `solve` (Muller on `det 𝒟`).

## 5. Scope / what is not done

- **Newberger closure** and the general `CoupledVDF` evaluator (the bicubic is fed in as a complex-analytic `f₀`).
- The outer `v⊥` integral is still Gauss–Kronrod. It too closes in finite form — the
  Bessel-weighted `∫ v⊥^d Jₙ₁Jₙ₂ dv⊥` per cell is the Bessel-product power series (Schläfli ₂F₃) in [`perp_analytic.jl`](../src/perp_analytic.jl) (`perp_pwpoly`) — which would make `χ_n` a pure finite double-cell sum (MPDES `intPar`×`intPer`), no quadrature at all.
  Not yet wired; `Jₙ′` must first be expanded as `(Jₙ₋₁−Jₙ₊₁)/2` to integer orders.
