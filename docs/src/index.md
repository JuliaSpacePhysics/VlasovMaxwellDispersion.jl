# VlasovMaxwellDispersion.jl

Linear hot-magnetized Vlasov–Maxwell dispersion solver, performant for arbitrary
velocity distribution functions (VDFs), analytic or numerical, with faster paths
for specialized cases (bi-Maxwellian, kappa, …).

## Quick start

```julia
using VlasovMaxwellDispersion

pl = NormalizedSpecies(-1.0, 1.0, Maxwellian(1.0))   # Ω̃, Π̃², VDF
k  = Wavenumber(0.0, 0.7)                            # k̃ = (k⊥, k∥)·c/Ω_ref

sol = solve(DispersionProblem(pl, 1.2 - 0.1im, k), Muller())        # seeded root
gsol = solve(GlobalDispersionProblem(pl, (0.5 - 0.6im, 2.5 + 0.1im), k))  # all roots in a box

ks = [Wavenumber(0.0, kz) for kz in 0.3:0.05:1.0]
ωs = solve(DispersionProblem(pl, 1.2 - 0.1im, ks)).omega            # track a k branch
```

The [Cattaert 2007 benchmark](cattaert.md) page works a full non-Maxwellian
benchmark end to end — seeded branch tracking, then a seedless global survey.
