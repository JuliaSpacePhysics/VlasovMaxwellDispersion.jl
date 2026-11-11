# VlasovMaxwellDispersion.jl

Linear hot-magnetized Vlasov–Maxwell dispersion solver. Finds roots of
`det 𝒟(ω,k) = 0` for gyrotropic distributions: analytic or tabulated,
non-relativistic or relativistic. Specialized distributions use closed-form
fast paths; arbitrary distributions use numerical integration.

```julia
using Pkg
Pkg.add(url = "https://github.com/JuliaPlasma/VlasovMaxwellDispersion.jl"); Pkg.add(["PlasmaPopulations", "Unitful"])
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
root = solve(DispersionProblem(plasma, (0.9 - 0.01im) * Ωp, k), Muller())
root.omega / Ωp, root.resid, root.retcode

# Find every root inside a complex-frequency box.
box = ((0.0 - 0.3im) * Ωp, (1.1 + 0.05im) * Ωp)
survey = solve(DispersionProblem(plasma, box, k), AAA())
survey.omega ./ Ωp
```

## Shared interface

Problems, targets, wavevectors, solutions, units, `track` and `dispersion_diagram` are
[PlasmaDispersionBase.jl](https://github.com/JuliaPlasma/PlasmaDispersionBase.jl#readme)'s, whose exported API this package re-exports; its README documents them.
Particles, `Population`, `Plasma` and the distribution models below come from
[PlasmaPopulations.jl](https://github.com/JuliaPlasma/PlasmaPopulations.jl#readme).

On top of those conventions (`e^{-iωt}`, so `Im ω > 0` is growth):

- Model speeds are m/s (bare or `Unitful`, e.g. `0.1u"c"`) in a physical `Plasma`, and `c` in a
  `NormalizedPopulation`; momenta of `SolverVDF`s (`CoupledVDF`, `SplineVDF`, …) are always in `mc`.
- Solvers run normalized, so numerical options (`link`'s `gate`, the `isgrowing` margin, `atol`) are in `Ω_ref`.

## Build plasma and distribution

To bypass physical units, give the normalized gyrofrequency `Ω̃`,
plasma-frequency squared `Π̃²`, and distribution, with speeds in `c`:

```julia
plasma = NormalizedPopulation(-1.0, 1.0, BiMaxwellian(; vth = 0.02))
```

Supported distribution:

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

## Algorithms

| target and `k` | algorithm |
|---|---|
| `Seed` at one `Wavevector` | `Muller()`, or a `Roots.jl` / `NonlinearSolve.jl` method once loaded |
| `Seed` along a path or a one-axis sweep | `Continuation(Muller())`: predicts along `k`, `base` corrects at each `k` |
| box (`Region`), any `k` | `AAA()` or `GRPF()` at each `k`; `link` joins the roots into branches |

```julia
# One branch, seeded at ks[10] and continued both directions.
solve(DispersionProblem(plasma, Seed(ω0, ks[10]), ks), Continuation(Muller()))

# Every branch in box over a sweep.
link(solve(DispersionProblem(plasma, (ω_lower_left, ω_upper_right), geometry), AAA()))
```

This solver's `solve` keywords, also storable on `DispersionProblem`:

| keyword | values |
|---|---|
| `mode` | `:det` (default); a polarization factor `:L`, `:R`, `:P`, `:O`, `:X`, see [Mode reduction](docs/src/reduction.md) |
| `backend` | `FixedNodeEval(; tol)` (default), fixed-node `χ` tables; `AdaptiveEval(; closure)`, per-ω adaptive quadrature with `HarmonicSum()` or `Newberger()` orbit integrals |

## Read results

Solutions and branches are described in the shared interface. Here `resid` is the scale-invariant
`σ_min/σ_max` of the full `𝒟` the solver evaluated; `residual(sol.prob, ω)` recomputes it on the
direct, unplanned path. See [return codes](docs/src/usage.md#return-codes-and-operation) for `retcode` meanings.

## More

- [Advanced usage](docs/src/usage.md): numerical controls, diagnostics
- [Case studies](docs/src/case-studies): executable published benchmarks
- [Mode reduction](docs/src/reduction.md) for parallel and perpendicular propagation: `:L`, `:R`, `:P`, `:O`, and `:X`
- [Physics derivation](docs/derivation.md)
- [Architecture](docs/architecture.md)
