# Coupled VDFs: choosing a representation

A gyrotropic `f₀(p⊥,p∥)` that does not factor as `f⊥(p⊥)·f∥(p∥)` — a spherical shell,
a bi-kappa, anything measured — is the expensive case. This page explains why, defines
the one quantity that decides the cost (the **separation rank**), compares the four
representations the package offers, and shows why the tabulated-spline route
([`GridVDF`](@ref)) silently returns wrong damping rates.

## Where the cost comes from

The non-relativistic magnetized susceptibility is

    χ(ω,k) ∝ Σₙ ∫dp⊥ ∫dp∥  Tₙ(k⊥p⊥/Ω) ⊗ ∇f₀(p⊥,p∥) / (p∥ − ζₙ),   ζₙ = (ω − nΩ)/k∥

Two structural facts drive everything:

1. **`ω` enters only through the pole `ζₙ`.** The Bessel bilinears `Tₙ`, the
   perpendicular quadrature nodes, the harmonic cap — none of them depend on `ω`. They
   should be built once per `k` and reused for every `ω`.
2. **For a coupled `f₀`, they cannot be.** The inner Cauchy transform
   `∫ ∂f₀(p⊥,p∥)/(p∥−ζₙ) dp∥` depends on `p⊥`, so the outer perpendicular integral wraps
   the `ω`-dependent inner one and the whole 2-D quadrature is redone at every `ω`.

That coupling is the entire problem. Measured for the proton shell at `θ=89.5°`,
`kλp=12.5`: 405 perpendicular nodes × (215 inner Landau evaluations + 75 harmonic
peel probes) = **118 000 `f₀` evaluations, ~42 ms — per `ω`**. An `AAA` survey asks for
~320 `ω` per `k`.

If instead `f₀ ≈ Σₛ ãₛ(p⊥)·bₛ(p∥)`, the perpendicular Bessel moments become
`ω`-independent tensors `P[n,s]` built once per `k`, and each `ω` costs only
`R·(2nmax+1)` **scalar** Cauchy transforms. The question is how big `R` has to be.

## Separation rank

View `f₀` as a bivariate kernel and sample it on its box. The **ε-rank** `R(ε)` is the
number of singular values of that sample matrix exceeding `ε·σ₁` — equivalently, the
smallest `R` for which some `Σₛ₌₁ᴿ ãₛ(p⊥)bₛ(p∥)` is within `ε` of `f₀`.

For analytic `f₀` the singular values decay geometrically, so `R(ε)` grows like
`log(1/ε)` — the rank is *small and barely sensitive to the tolerance*:

| `f₀` | `R(1e-4)` | `R(1e-6)` | `R(1e-8)` | `R(1e-10)` |
|---|---|---|---|---|
| bi-Maxwellian (separable) | 1 | 1 | 1 | 1 |
| skewed Gaussian `exp[−(p∥²+p⊥²+0.6p∥p⊥)]` | 4 | 6 | 8 | 10 |
| proton shell `exp[−(√(p⊥²+p∥²)−v_d)²/c_p²]` | 7 | 9 | 11 | 13 |
| bi-kappa, `κ=3` | 4 | 7 | 11 | 14 |

So the hard cases are rank ~10. That is the whole reason this works.

**Rank is bought with conditioning.** A rank-`R` cross (skeleton) factorization inverts
the `R×R` pivot matrix `f₀[vₚ,uₚ]`, whose condition number is `≈1/rtol`. Pushing `rtol`
down raises `R` slowly but degrades the perpendicular factors `ãₛ` proportionally, so the
achievable accuracy in `χ` floors at roughly `100·rtol`. In practice `rtol=1e-8` (the
default) gives `χ` to ~1e-5–1e-9 depending on the VDF; `rtol=1e-10` buys another decade.
This also means **you must never run an adaptive quadrature to a tolerance tighter than
`rtol` on the `ãₛ`** — it will chase round-off and subdivide without bound.

## The four representations

### Separable families — [`Maxwellian`](@ref), [`ProductBiKappa`](@ref), [`SeparableVDF`](@ref)

The parallel Cauchy moments close analytically for Gaussian and kappa factors; arbitrary
[`SeparableVDF`](@ref) factors use a Landau quadrature. Perpendicular Bessel moments are
cached once per fixed `k`.

- **Pro:** correct within the configured truncation and quadrature tolerance at any `Im ω`
  where the factors are analytic.
- **Con:** only for distributions that factor as `f⊥(p⊥)f∥(p∥)`.

### [`CoupledVDF`](@ref) — adaptive 2-D quadrature

- **Pro:** exact for any analytic `f₀`. The Landau residue is taken from `f₀` itself, so
  the analytic continuation to damped `ω` is the true one. This is the **reference**.
- **Con:** ~1e5 `f₀`-evaluations per `ω`; 2.6–42 ms for Case 5. Unusable for a survey.

### [`GridVDF`](@ref) — tensor spline of tabulated data

- **Pro:** the only option when `f₀` *is* data (spacecraft, PIC output). `f₀` becomes
  piecewise-polynomial, so the parallel Hilbert transform closes per cell.
- **Con:** slow in practice (891 ms/eval on the Case-5 shell — the harmonic loop sits
  outside the perpendicular quadrature), and **it cannot be analytically continued**. See
  below.

### [`LowRankVDF`](@ref) — adaptive-cross skeleton

`f₀ ≈ Σₛ ãₛ(p⊥)·bₛ(p∥)` with `bₛ(u) = f₀(vₛ,u)` and `ãₛ(v) = Σᵣ f₀(v,uᵣ)·M[r,s]`.

- **Pro:** `O(R·nmax)` per `ω`; 117–313× faster than `CoupledVDF` on Case 5 (22–134 µs).
  Construction ~8 ms, per-`k` plan ~2 ms. Crucially, the parallel factors are **literal
  slices of the true `f₀`**, so the Landau residue is exact and damped roots are as
  accurate as the reference path.
- **Con:** approximate (`isexact` is `false`); accuracy is capped by `rtol` and its
  conditioning; needs a bounded box, so heavy tails must be truncated.

A cross approximation is used rather than an SVD precisely because an SVD's factors are
numerical vectors that exist only on the sample grid — they have no continuation off the
real axis. The cross approximates *the coupling between the two coordinates*, not the
analytic structure of either one.

## Why the spline fails for damped roots

Split the Landau-continued transform into the two pieces the algorithm actually computes:

    ∫ g(u)/(u−ζ) du │continued from Im ζ>0  =  ∫ g(u)/(u−ζ) du  +  σ·2πi·g(ζ)
                                                └── on the real axis ──┘   └─ at a COMPLEX point ─┘

Replace the true `g` by a surrogate `ĝ` with `‖g−ĝ‖∞ = ε` on the real axis.

- The **first term is stable**: its error is bounded by `ε·∫|du|/|u−ζ| = ε·O(log)`. No
  amplification. A spline is fine here.
- The **second term is not**. It needs `g` at `ζ`, a distance `d = |Im ω| / |k∥|` *off*
  the real axis, where the surrogate was never fitted. A degree-`p` polynomial on a cell
  of width `h` amplifies its fit error by roughly `(2d/h)^p` when continued a distance
  `d`, while the true `g` follows its own analytic structure (for a Gaussian tail,
  `exp(d²/c²)`). The two have nothing to do with each other.

The controlling ratio is `d/h = |Im ω| / (|k∥|·h∥)`, cells deep. For Case 5 with a
161×321 grid: `h∥ ≈ 8.8e-5`, `k∥ = 32.3`, and at `Im ω = −0.06`, `d = 1.9e-3` — **21 cells
below the axis**, extrapolated with a cubic. Measured relative error of `GridVDF`'s `χ`
against the exact path (with a tight `rtol=1e-7` fit, so this is *not* fit error):

| `Im ω` | `+0.02` | `0` | `−0.005` | `−0.02` | `−0.06` |
|---|---|---|---|---|---|
| `GridVDF` rel. error | 1.2e-3 | 3.0e-3 | 1.8e-4 | 6.4e-3 | **0.90** |
| `LowRankVDF` rel. error | — | — | 3e-9 | 2e-9 | 4e-8 |

Two things worth internalizing:

1. **Refining the grid makes it worse, not better.** Refining shrinks `h`, which *raises*
   `d/h`. There is no grid at which a compactly-supported real-axis surrogate learns the
   continuation. This is structural, not a tuning failure.
2. It is silent. `GridVDF` returns a finite, plausible `χ`; only the imaginary part is wrong.

**Practical rule:** `GridVDF` is for real or very weakly damped `ω` (`|Im ω| ≲ |k∥|·h∥`).
If `f₀` is tabulated *and* damping matters, fit an analytic model first (a sum of
Gaussians or kappas, whose continuations are known in closed form) and feed that to
`LowRankVDF` or a sum of [`SeparableVDF`](@ref)s.

## The damping wall (physics, not numerics)

Even with an exact continuation there is a floor on how deep a root can be sought. The
Landau continuation of a Gaussian-tailed `f₀` grows like `exp[(Im ω/(k∥·c))²]`. Measured
for the Case-5 shell at `kλp=12.5`:

| `Im ω` | `−0.005` | `−0.06` | `−0.12` | `−0.2` | `−0.35` |
|---|---|---|---|---|---|
| `‖χ‖` | 8.2e2 | 2.8e3 | 4.2e5 | 5.8e10 | 1.2e27 |

Past `Im ω ≈ −0.4` even `CoupledVDF` overflows to `NaN`. The usable box depth is
`|Im ω| ≲ few × k∥·v_th`; at `θ=89.5°` that is `k∥c_p ≈ 0.047`, which is exactly where the
Case-5 region's floor (`−0.06`) sits. Propagating closer to perpendicular buys harmonics
but *costs* accessible damping.

For kappas the limit is sharper and algebraic: the parallel slice
`(1 + p∥²/a∥ + p⊥²/a⊥)^{-(κ+1)}` has branch points at `|Im p∥| = √(a∥(1+p⊥²/a⊥)) ≥ √a∥`,
so no continuation exists at all beyond `|Im ω| < |k∥|·√a∥`.

## Choosing

| situation | use |
|---|---|
| `f₀` is in a closed-form family | that family — `Maxwellian`, `ProductBiKappa`, `BiKappa` |
| `f₀` factors | [`SeparableVDF`](@ref) |
| coupled analytic `f₀`, want speed (surveys) | [`LowRankVDF`](@ref) |
| coupled analytic `f₀`, want a certified reference | [`CoupledVDF`](@ref) |
| tabulated `f₀`, real / weakly damped `ω` | [`GridVDF`](@ref) |
| tabulated `f₀`, damped roots matter | fit an analytic model, then `LowRankVDF` |

Cross-check by construction: `LowRankVDF` of a separable `f₀` returns rank 1 and
reproduces the bi-Maxwellian closed form to ~1e-12, including at `Im ω = −2`.
