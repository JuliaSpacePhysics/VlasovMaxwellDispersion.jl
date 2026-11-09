# VlasovMaxwellDispersion.jl

Linear hot-magnetized Vlasov–Maxwell dispersion solver. Finds roots of
`det 𝒟(ω,k) = 0` for gyrotropic distributions: analytic or tabulated,
non-relativistic or relativistic. Specialized distributions use closed-form
fast paths; arbitrary distributions use numerical integration.

```julia
using Pkg
Pkg.add(url = "https://github.com/JuliaSpacePhysics/VlasovMaxwellDispersion.jl"); Pkg.add(["PlasmaPopulations", "Unitful"])
```

## Quick start

```julia
using VlasovMaxwellDispersion, Unitful
using PlasmaPopulations: gyrofrequency, inertial_length

p, e = Proton(), Electron()
n = 5u"cm^-3"
plasma = Plasma(
    Population(p, BiMaxwellian; T_perp = 100u"eV", T_para = 20u"eV", n),
    Population(e, BiMaxwellian; T = 50u"eV", n);
    B0 = 5u"nT",
)

Ωp = gyrofrequency(p, plasma.B0)        # rad/s
k = Wavevector(0.0, 0.5 / inertial_length(p, n))  # k∥ dᵢ = 0.5, in m⁻¹

# Refine one root from a nearby complex frequency.
root = solve(DispersionProblem(plasma, (0.9 - 0.01im) * Ωp, k))
root.omega / Ωp, root.resid, root.retcode

# Find every root inside a complex-frequency box.
box = ((0.0 - 0.3im) * Ωp, (1.1 + 0.05im) * Ωp)
survey = solve(DispersionProblem(plasma, box, k))
survey[1].omega / Ωp
```

## Conventions

Units follow the plasma you pass. `Ω_ref` is a reference gyrofrequency, the first population's by default.

| quantity | physical `Plasma` | `NormalizedPlasma`, `NormalizedPopulation` | note |
|---|---|---|---|
| `ω`, boxes, seeds | rad/s | `Ω_ref` | complex; `e^{-iωt}`, so `Im ω > 0` means growth |
| `k` | m⁻¹ | `Ω_ref/c` | `Wavevector(k_perp, k_para)`, with `B₀ ∥ ẑ` |
| model speeds | m/s, bare or `Unitful` (`0.1u"c"`) | `c` | |
| momenta of `SolverVDF`s (`CoupledVDF`, `SplineVDF`, …) | `mc` | `mc` | |

Solvers always run normalized, so their numerical options (`linking.gate`, the `isgrowing` margin, `atol`) are in `Ω_ref`.

## Build plasma and distribution

Particles (`Particle`, `Electron`, `Positron`, `Proton`), `Population`, `Plasma`, the
distribution models below and `Wavevector`, `AngleSweep`, `CartesianSweep` come from
[PlasmaPopulations.jl](https://github.com/JuliaPlasma/PlasmaPopulations.jl#readme) and are re-exported.
See its README for constructors and conventions.

`DispersionProblem(plasma, target, k; ref)` keeps SI inputs and results and solves internally
with `Ω_ref` from `ref`, another particle's gyrofrequency or a frequency [rad/s]. 

To work in normalized units, pass `NormalizedPlasma(plasma; ref)`:

```julia
DispersionProblem(NormalizedPlasma(plasma; ref = Electron()), box_in_Ωe, k_in_Ωe_over_c)
```

`scales(plasma; ref)` gives that normalization's dimensionless scales: `vA`,
`Omega_ref`, and per population `wc`, `wp`, inertial length `d`, thermal speeds and gyroradius `rho`.

To bypass physical units, give the normalized gyrofrequency `Ω̃`,
plasma-frequency squared `Π̃²`, and distribution, with speeds in `c`:

```julia
plasma = NormalizedPopulation(-1.0, 1.0, BiMaxwellian(; vth = 0.02))
```

Distribution choices:

| constructor | use |
|---|---|
| `Cold()` | cold limit (closed Stix form) |
| `BiMaxwellian(; vth_perp, vth_para, vd_para)` | drifting bi-Maxwellian |
| `GaussianRing(; vr_perp, w_perp, vth_para, vd_para)` | shifted-Gaussian ring (beam with `vd_para`), accurate for `k⊥ vr_perp/Ω ≲ 10` |
| `BiKappa(; vth_perp, vth_para, kappa)` | coupled bi-kappa, `κ > 3/2` |
| `ProductBiKappa(; vth_perp, vth_para, kappa_perp, kappa_para)` | separable product bi-kappa, no loss cone |
| `GyroRing(vth_perp, vr) ⊗ Gaussian(vth_para, vd)` | gyro-averaged perpendicular-drifting Maxwellian (`I₀` ring) |
| `MaxwellJuttner(; theta)` | isotropic relativistic distribution, `theta = T/(mc²)` |
| `SeparableVDF(fperp, fpara; perp, para)` | analytic `f₀ = f⊥ f∥` |
| `CoupledVDF(f0; perp, para)` | arbitrary analytic `f₀(p⊥,p∥)`, e.g. `(q, u) -> pdf(model, q, u)` for a `model` with speeds in `c` |
| `LowRankVDF(f0; perp, para)` | faster surrogate for a coupled `f₀` |
| `GridVDF(v_perp, v_para, f)` | tabulated or simulation data, fitted by `SplineVDF`'s defaults |
| `SplineVDF(vperp, vpara, F; method, regime)` | the fit itself, momenta in `mc`; relativistic data |
| `ReducedVDF(fpara; para)` | 1-D parallel `f∥` for field-aligned electrostatic roots (`k⊥ = 0`) |

Kappa models do not support a drift `vd_para`. Analytic functions must accept complex
parallel arguments for Landau continuation. Use `GridVDF` for real-axis-only
data. Prefer `CoupledVDF` over `LowRankVDF` for strongly damped roots.

See [Advanced usage](docs/src/usage.md) for constructor details, integration
bounds, custom gradients, relativistic coordinates.

## Choose wavevectors and roots

| `k` argument | result |
|---|---|
| `Wavevector(k_perp, k_para)` | one wavevector |
| vector of `Wavevector`s | one seeded branch continued along path |
| `AngleSweep(; k, theta)` | polar sweep |
| `CartesianSweep(; k_perp = 0, k_para)` | Cartesian sweep |

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

A seeded solve returns `DispersionSolution`, in the units of the problem:

| field | meaning |
|---|---|
| `.omega` | root, or roots along path; failed points are `NaN` |
| `.resid` | scale-invariant `σ_min/σ_max` of full `𝒟` |
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
