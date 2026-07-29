# Advanced usage

See [README](https://github.com/JuliaSpacePhysics/VlasovMaxwellDispersion.jl#readme) for conventions, plasma construction, distribution selection and result shapes. 
This page covers controls beyond that default workflow.

## Constructor details

For `SeparableVDF`, `CoupledVDF`, and `LowRankVDF`, integration bounds may be
explicit tuples or bare half-widths:

```julia
CoupledVDF(f0; perp = (0.0, 6.0), para = (-8.0, 4.0))
CoupledVDF(f0; perp = 6.0, para = 6.0)  # (0, 6) × (-6, 6)
```

Use explicit shifted bounds for drifting or ring-shifted functions. Gradients
are differentiated automatically; override them with `dgrad`, or `dfperp` and
`dfpara`.

General representations:

```julia
# Separable analytic function.
κ, vth = 6.0, 0.05
fκ(x) = (1 + x^2 / (κ * vth^2))^(-κ)
separable = SeparableVDF(fκ, fκ; perp = 12vth, para = 12vth)

# Arbitrary coupled function and faster real-axis surrogate.
f0(q, u) = exp(-(u^2 + q^2 + 0.6u * q))
coupled = CoupledVDF(f0; perp = 6.0, para = 6.0)
lowrank = LowRankVDF(f0; perp = 6.0, para = 6.0, rtol = 1e-8)

# Tabulated data.
vperp, vpara = range(0, 6, 61), range(-6, 6, 81)
F = [exp(-(v^2 + u^2)) for v in vperp, u in vpara]
grid = GridVDF(vperp, vpara, F)
```

Closed-form one-dimensional factors compose with `⊗`:

```julia
vdf = Gaussian(0.02) ⊗ Kappa(0.02, 4)
```

This combines a Maxwellian perpendicular core with a parallel kappa tail.
Available factors include `Gaussian(vth, vd)`, `GyroRing(vth, vr)`, and
`Kappa(θ, κ)`.

Relativistic energy coordinates take `f0(γ,p∥)` and
`dgrad(γ,p∥) = (∂γ f, ∂∥ f)`:

```julia
CoupledVDF(
    f0; regime = Relativistic(), perp, para,
    dgrad, coords = :energy
)
```

Representation-specific behavior lives in:

- [Choosing a coupled representation](representations/coupled-vdf.md)
- [GridVDF projection](representations/grid-vdf.md)
- [Low-rank susceptibility](representations/lowrank-susceptibility.md)

## Sweep geometry

Each `AngleSweep` or `CartesianSweep` axis may be a fixed scalar, explicit
sample vector, or `(lo, hi)` tuple. Tuples generate 61 samples.

```julia
ks = [Wavenumber(0.0, kz) for kz in 0.3:0.05:1.0]
angle = AngleSweep(k = (0.01, 0.6), theta = deg2rad(45))
cartesian = CartesianSweep(kz = 0.01:0.008:1.0)
surface = CartesianSweep(kx = (0.0, 1.0), kz = (0.0, 1.0))
```

Surveys accept two swept axes. Seeded continuation accepts one.

## Algorithms

Defaults follow target and geometry:

| target | default | main controls |
|---|---|---|
| seed at fixed `k` | `Muller()` | `atol`, `maxiter` |
| seed along path | `Continuation()` | `base`, `order`, `reltol`, `abstol`, `maxsubdiv` |
| frequency box | `AAA()` | `n = (nRe,nIm)`, `tol`, `max_degree` |

Override default by passing an algorithm:

```julia
solve(problem, AAA(n = (80, 80)))
```

Survey-only keywords control polishing and branch linking:

```julia
solve(problem; refine = Muller(), linking = (;))
```

`mode = :det` is valid for every geometry. Cheaper polarization factors have
strict symmetry and propagation constraints; see [Mode reduction](reduction.md)
before selecting `:L`, `:R`, `:P`, `:O`, or `:X`.

## Susceptibility backends

`FixedNodeEval(; order, nprobe)` is default. It plans quadrature nodes at each
wavevector and reuses frequency-independent work. Use `AdaptiveEval()` when
fixed nodes under-resolve a distribution; it repeats adaptive quadrature at
every frequency and is slower.

```julia
problem = DispersionProblem(
    plasma,
    target,
    k;
    backend = AdaptiveEval(),
)
```

## Return codes and operation

| code | meaning |
|---|---|
| `Success` | every requested root converged |
| `Partial` | branch tracked, but some wavevectors contain `NaN` |
| `Saturated` | survey fit stopped before tolerance and may have missed roots |
| `MaxIters` | root polisher reached iteration limit |
| `Failure` | no root found |

For `Saturated`, shrink frequency box or increase `AAA(n = ...)`.
`successful_retcode(solution)` tests for full convergence.

Surveys are threaded over wavevector points; start Julia with `-t auto`.
`Region` boundaries are soft, so branches may continue slightly past box edges.

## Diagnostic API

For custom root finding, diagnostics, or determinant maps:

```julia
𝒟(plasma, ω, k)
residual(plasma, ω, k)

# Callable scalar function with frequency-independent work planned once.
D = DispersionFunction(plasma, k)

# Per-species susceptibility.
species = NormalizedPlasma(plasma).species[1]
contribution(species, ω, k)

# Bare distribution at Ω̃ = Π̃² = 1.
contribution(vdf, ω, k)
```
