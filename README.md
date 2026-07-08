# VlasovMaxwellDispersion.jl

Linear hot-magnetized Vlasov–Maxwell dispersion solver. Performant for arbitrary velocity distribution functions (VDFs), analytic or numerical. Even faster for specialized cases (e.g., bi-Maxwellian, kappa, etc.).

Full derivation: [`derivation.md`](docs/derivation.md).
Method and internals: [`architecture.md`](docs/architecture.md).

## Quick start

```julia
using Pkg; Pkg.develop(path="."); Pkg.instantiate()
using VlasovMaxwellDispersion

pl = NormalizedSpecies(-1.0, 1.0, Maxwellian(1.0))   # Ω̃, Π̃², VDF
k  = Wavenumber(0.0, 0.7)                           # k̃ = (k⊥, k∥)·c/Ω_ref

alg = Muller()
sol = solve(DispersionProblem(pl, 1.2 - 0.1im, k), alg)             # seeded root (Langmuir+Landau)
ω   = sol.omega

gsol = solve(GlobalDispersionProblem(pl, (0.5 - 0.6im, 2.5 + 0.1im), k))  # all roots in a box
roots = gsol.roots

ks   = [Wavenumber(0.0, kz) for kz in 0.3:0.05:1.0]
ωs   = solve(DispersionProblem(pl, 1.2 - 0.1im, ks)).omega          # track k branch
```

## Distributions

Every VDF below plugs into the same `solve` workflow — only the constructor
changes. `Ω̃`, `Π̃²` are the per-species inputs to `NormalizedSpecies`; speeds are
`v/c`.

```julia
# electron–proton, all referenced to first species' frequency |Ω_e|
mp_me = 1836.15

# bi-Maxwellian: temperature anisotropy + parallel drift (whistler / firehose)
pl = NormalizedSpecies(-1.0, 4.0, Maxwellian(vth_para=0.02, vth_perp=0.04, vd=0.0))
ω  = solve(DispersionProblem(pl, 0.9 - 0.05im, Wavenumber(0.0, 0.5))).omega

# Ring / shell: perp ring speed vr — I₀ gyro-ring, or a literal shifted Gaussian
ring  = Maxwellian(vth_para=0.1, vth_perp=0.1, vr=0.3)
shell = GaussianRing(vth_para=0.1, vth_perp=0.1, vr=0.3)

# Field-aligned electrostatic (k⊥=0): bump-on-tail / two-stream / Landau
f(u) = 0.94exp(-u^2)/√π + 0.06exp(-((u-4)/0.5)^2)/(0.5√π)
pl = NormalizedSpecies(-1.0, 1.0, ReducedVDF(f; para=(-12.0, 14.0)))
ω  = solve(DispersionProblem(pl, 1.0 + 0.05im, Wavenumber(0.0, 0.25))).omega  # Im ω>0 ⇒ growth

# Arbitrary f₀(p⊥,p∥)
g(q, u) = exp(-(u^2 + q^2 + 0.6u*q))
pl = NormalizedSpecies(-1.0, 1.0, CoupledVDF(g; perp=6.0, para=(-6.0, 6.0)))

# Arbitrary separable f₀ = f⊥(p⊥)·f∥(p∥) — e.g. a product-bi-kappa
κ, vth = 6.0, 0.05
fκ(x) = (1 + x^2/(κ*vth^2))^(-κ)
pl = NormalizedSpecies(1.0, 1.0, SeparableVDF(fκ, fκ; perp=12vth, para=(-12vth, 12vth)))

# Tabulated / numerical VDF on a (v⊥, v∥) grid → positivity-preserving NNLS B-spline
vperp, vpar = range(0, 6, 61), range(-6, 6, 81)
F = [exp(-(v^2 + u^2)) for v in vperp, u in vpar]   # F[perp, para]
pl = NormalizedSpecies(-1.0, 0.5, GridVDF(vperp, vpar, F))

# Relativistic isotropic Maxwell–Jüttner (μ = mc²/T)
pl = NormalizedSpecies(-1.0, 0.5, MaxwellJuttner(mu=40.0))
```

### From physical units

Build species physically (SI or `Unitful` via the extension) and normalize to reference.
`NormalizedSpecies(s::Species, B0, ref)` maps `(B0, n, q, m) →(Ω̃, Π̃²)`; `k` and `ω` stay dimensionless (`·c/Ω_ref`, `·/Ω_ref`).

```julia
using VlasovMaxwellDispersion.PlasmaBase   # Electron, Proton, Species

B0 = 5.0e-9                                 # Tesla
e  = Species(Electron(), Maxwellian(0.02); n=5.0e6)
p  = Species(Proton(),   Maxwellian(0.02/√1836); n=5.0e6)
plasma = (NormalizedSpecies(e, B0, Proton()), NormalizedSpecies(p, B0, Proton()))
```

## Capabilities

| VDF / mode | constructor | method |
|---|---|---|
| Arbitrary `f₀(p⊥,p∥)` | `CoupledVDF` | nested 2-D orbit quadrature |
| Arbitrary separable `f⊥(p⊥)·f∥(p∥)` | `SeparableVDF` | Hilbert × Bessel-moment quadrature |
| Grid / numerical | `GridVDF` | NNLS B-spline → piecewise-poly `H∥`,`P⊥` |
| 1-D parallel (electrostatic, `k⊥=0`) | `ReducedVDF` | parallel Hilbert transform |

Specialized VDFs include: Cold fluid `ColdVDF`, (bi-)Maxwellian, drifting, ring `Maxwellian` / `GaussianRing`, Maxwell–Jüttner (relativistic) `MaxwellJuttner`.

Solvers follow the `CommonSolve.solve(problem, algorithm)` interface: a
`DispersionProblem`/`GlobalDispersionProblem` (seeded / seedless survey) solved by
`Muller` (default) / `GRPF` (`RootsAndPoles.jl`) / `ArcLength`. `solve` returns a
`DispersionSolution` (`.omega`, `.retcode`).

Two closures for the orbit integral are available (`derivation.md` §3), passed as
the `closure=` keyword of `solve`:

- `HarmonicSum()` (default): truncate the harmonic sum at `nmax ≈ k⊥ρ`. Handles damping via the Landau contour.
- `Newberger()`: Qin closed-orbit `T(a,z)` in complex-order Bessel `J_{±a}`, no truncation, damping via residue extraction.
