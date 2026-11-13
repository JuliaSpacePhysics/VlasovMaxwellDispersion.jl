# Advanced usage

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
grid = SplineVDF(vperp, vpara, F)
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
- [SplineVDF projection](representations/grid-vdf.md)
- [Low-rank susceptibility](representations/lowrank-susceptibility.md)

## Algorithms

Surveys accept sweeps with two swept axes; `Continuation` accepts one.

| algorithm | main controls |
|---|---|
| `Muller()` | `atol`, `maxiter` |
| `Continuation(base)` | `order`, `reltol`, `abstol`, `maxsubdiv` |
| `AAA()`, `GRPF()` | `n = (nRe,nIm)`, `tol`, `max_degree`, `refine` |

```julia
solve(problem, AAA(n = (80, 80)))
```

`AAA` and `GRPF` polish their roots with `refine` (`Muller()` default, `nothing` keeps raw roots);
`link` joins a survey's roots into branches:

```julia
survey = solve(problem, AAA(refine = Muller(atol = 1e-12)))
link(survey; gate = 0.05)
```

`mode = :det` is valid for every geometry. Cheaper polarization factors have
strict symmetry and propagation constraints; see [Mode reduction](reduction.md)
before selecting `:L`, `:R`, `:P`, `:O`, or `:X`.

## Susceptibility backends

`FixedNodeEval(; order, nprobe, tol)` is default. It tabulates `∇f₀` on fixed nodes to
interpolation tolerance `tol` and reuses frequency-independent work. Use `AdaptiveEval()` when
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

Codes are `ReturnCode`'s (see its docstring). Here `Saturated` means the `AAA` fit or `GRPF` mesh
stopped short of its tolerance and may have missed roots: shrink the frequency box or increase
`AAA(n = ...)`. `Partial` means `Continuation` left some wavevectors `NaN`.

Surveys are threaded over wavevector points; start Julia with `-t auto`.

## Diagnostic API

For custom root finding, diagnostics, or determinant maps:

```julia
dispersion_tensor(plasma, ω, k)   # 𝒟, also `VlasovMaxwellDispersion.𝒟`
residual(plasma, ω, k)            # direct path, independent of the plan behind `sol.resid`
electrostatic_det(plasma, ω, k)   # k̂·ε·k̂, zero at electrostatic modes

# Callable scalar function with frequency-independent work planned once; normalized plasma only.
D = VlasovMaxwellDispersion.DispersionFunction(NormalizedPlasma(plasma), k)

# Per-population susceptibility.
pop = NormalizedPlasma(plasma).populations[1]
contribution(pop, ω, k)

# Bare distribution at Ω̃ = Π̃² = 1.
contribution(vdf, ω, k)
```

These convert each population's model to this solver's form on every call. For repeated
evaluation, convert once with `prepare(plasma)`; for a `GridVDF`, conversion is a spline fit.
