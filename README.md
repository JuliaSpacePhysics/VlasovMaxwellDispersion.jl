# VlasovMaxwellDispersion.jl

Linear hot-magnetized Vlasov–Maxwell dispersion solver. Finds roots of
`det 𝒟(ω,k) = 0` for gyrotropic distributions: analytic or tabulated,
non-relativistic or relativistic. Specialized distributions use closed-form
fast paths; arbitrary distributions use numerical integration.

```julia
using Pkg
Pkg.add(url = "https://github.com/JuliaSpacePhysics/VlasovMaxwellDispersion.jl"); Pkg.add("Unitful")
```

## Quick start

```julia
using VlasovMaxwellDispersion, Unitful

plasma = Plasma(
    Species(Proton(); n = 5u"cm^-3", T = (100u"eV", 20u"eV")),
    Species(Electron(); n = 5u"cm^-3", T = 50u"eV");
    B0 = 5u"nT",
)

s = scales(plasma, 1)
k = Wavenumber(0.0, 0.5 / s.d)  # k∥ dᵢ = 0.5

# Refine one root from a nearby complex frequency.
root = solve(DispersionProblem(plasma, 0.9 - 0.01im, k))
root.omega, root.resid, root.retcode

# Find every root inside a complex-frequency box.
box = (0.0 - 0.3im, 1.1 + 0.05im)
survey = solve(DispersionProblem(plasma, box, k))
survey[1].omega
```

Frequencies and wavevectors are normalized to a reference gyrofrequency
`Ω_ref`, the first species' by default. Above, `ω` is in proton gyrofrequencies
and `k` is in `Ω_ref/c`.

## Conventions

| quantity | solver unit | note |
|---|---|---|
| `ω`, boxes, seeds | `Ω_ref` | complex; `e^{-iωt}`, so `Im ω > 0` means growth |
| `k` | `Ω_ref/c` | `Wavenumber(kperp, kz)`, with `B₀ ∥ ẑ` |
| `vth`, `vd`, `vr` | `c` | |
| physical `n`, `T`, `B0` | m⁻³, eV, tesla | bare numbers or `Unitful` quantities |

Use `Plasma(...; ref = particle)` to change `Ω_ref`. `scales(plasma, i)` gives
dimensionless scales for species `i`, including inertial length `d`, gyroradius
`rho`, thermal speeds, gyrofrequency ratio `wc`, and plasma-frequency ratio
`wp`.

Scalar thermal parameters are isotropic. A pair means `(⊥, ∥)`:
`T = (100u"eV", 20u"eV")` or
`T = (perp = 100u"eV", para = 20u"eV")`.

## Build plasma and distribution

Physical input:

```text
Species(particle, vdf = nothing; n, T | beta | vth)
Plasma(species...; B0, ref = first species' particle)
```

`particle` may be `Electron()`, `Proton()`, `Particle(; z, A)`, or
`Particle(q, m)`. Omit `vdf` for a bi-Maxwellian built from `T`, `beta`, or
`vth`; omit thermal input too for a cold species.

Choose another distribution with a function of species scales:

```julia
Species(
    Proton(),
    sc -> ProductBiKappa(sc; kappa = (200.0, 3.0));
    n = 2.4u"cm^-3",
    T = 2555u"eV",
)
```

Or bypass physical units by supplying normalized gyrofrequency `Ω̃`,
plasma-frequency squared `Π̃²`, and distribution:

```julia
plasma = NormalizedSpecies(-1.0, 1.0, Maxwellian(0.02))
```

Main distribution choices:

| constructor | use |
|---|---|
| `ColdVDF()` | cold limit |
| `Maxwellian(; vth, vd, vr)` | drifting bi-Maxwellian or gyro-ring |
| `GaussianRing(; vth, vd, vr)` | shifted-Gaussian ring; accurate for `k⊥vr/Ω ≲ 10` |
| `BiKappa(; vth, kappa)` | coupled bi-kappa, `κ > 3/2` |
| `ProductBiKappa(; vth, kappa)` | separable product bi-kappa |
| `MaxwellJuttner(; mu)` | isotropic relativistic distribution, `mu = mc²/T` |
| `SeparableVDF(fperp, fpara; perp, para)` | analytic `f₀ = f⊥ f∥` |
| `CoupledVDF(f0; perp, para)` | arbitrary analytic `f₀(p⊥,p∥)` |
| `LowRankVDF(f0; perp, para)` | faster surrogate for a coupled `f₀` |
| `GridVDF(vperp, vpara, F)` | tabulated or simulation data |

Analytic functions must accept complex parallel arguments for Landau
continuation. Use `GridVDF` for real-axis-only data. Prefer `CoupledVDF` over
`LowRankVDF` for strongly damped roots.

See [Advanced usage](docs/src/usage.md) for constructor details, integration
bounds, custom gradients, relativistic coordinates.

## Choose wavevectors and roots

| `k` argument | result |
|---|---|
| `Wavenumber(kperp, kz)` | one wavevector |
| vector of `Wavenumber`s | one seeded branch continued along path |
| `AngleSweep(; k, theta)` | polar sweep |
| `CartesianSweep(; kx = 0, kz)` | Cartesian sweep |

Target determines root search:

```julia
# One branch, seeded at ks[10] and continued both directions.
solve(DispersionProblem(plasma, Seed(ω0, ks[10]), ks))

# Every branch in box over a sweep.
solve(DispersionProblem(plasma, (ω_lower_left, ω_upper_right), geometry))
```

A number or `Seed` refines one branch. A complex box or `Region` surveys every
branch found inside it. Default algorithms are `Muller` for one point,
`Continuation` for a path, and `AAA` for a survey.

## Read results

A seeded solve returns `DispersionSolution`:

| field | meaning |
|---|---|
| `.omega` | root, or roots along path; failed points are `NaN` |
| `.resid` | scale-invariant `σ_min/σ_max` residual |
| `.retcode` | convergence status; see [return codes](docs/src/usage.md#return-codes-and-operation) |
| `.stats` | function evaluations and elapsed time |

A survey returns `SurveySolution`, a vector of branches:

```julia
branch = survey[1]
branch.omega, branch.k
```

## Plotting

```julia
using CairoMakie

fig, (axr, axi), plots = dispersion_diagram(survey; title = "…")
dispersion_diagram!((axr, axi), filter(isgrowing, survey))
```

## More

- [Advanced usage](docs/src/usage.md): numerical controls, plotting, diagnostics
- [Case studies](docs/src/case-studies): executable published benchmarks
- [Mode reduction](docs/src/reduction.md) for parallel and perpendicular propagation: `:L`, `:R`, `:P`, `:O`, and `:X`
- [Physics derivation](docs/derivation.md)
- [Architecture](docs/architecture.md)
