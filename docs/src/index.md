# VlasovMaxwellDispersion.jl

Linear hot-magnetized Vlasov–Maxwell dispersion solver, performant for arbitrary
velocity distribution functions (VDFs), analytic or numerical, with faster paths
for specialized cases (bi-Maxwellian, kappa, …).

## Quick start

```julia
using VlasovMaxwellDispersion, Unitful

plasma = Plasma(
    Species(Proton(); n = 5u"cm^-3", T = (100u"eV", 20u"eV")),  # T = (T⟂, T∥)
    Species(Electron(); n = 5u"cm^-3", T = 50u"eV"),
    B0 = 5u"nT",                                  # ω, k in Ω_ref = the proton gyrofreq
)

s = scales(plasma, 1)   # dimensionless d, ρ, vth_perp, vth_para, vA, … for species 1
k = Wavenumber(0.0, 0.5 / s.d)           # k∥·dᵢ = 0.5

sol = solve(DispersionProblem(plasma, 0.9 - 0.01im, k))      # seeded root
gsol = solve(DispersionProblem(plasma, (0.5 - 0.6im, 2.5 + 0.1im), k))  # all roots in a box
```

The [Cattaert 2007 benchmark](case-studies/cattaert.md) page works a full non-Maxwellian
benchmark end to end — seeded branch tracking, then a seedless global survey.
