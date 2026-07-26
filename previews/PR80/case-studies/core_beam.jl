# # Ion cyclotron waves from a proton core-beam
# 
# Counter-propagating ion cyclotron waves excited by a proton core-beam
# distribution. Both proton populations are drifting bi-Maxwellians with
# `T⊥ > T∥`; the combined anisotropy and beam drift destabilize forward- and
# backward-propagating ICWs with unequal growth rates.

# Reference: Case 6 in Guo et al. (2026, arXiv:2606.14439)

using VlasovMaxwellDispersion
using VlasovMaxwellDispersion: isgrowing
using Unitful

# ## Plasma setup
#
# Drift speeds are in `c`.

B0 = 7.5e-7u"T"
plasma = Plasma(
    Species(Proton(), sc -> Maxwellian(sc; vd=-2.87e-4); n=2.53e9u"m^-3", T=(20, 100) .* u"eV"),
    Species(Proton(), sc -> Maxwellian(sc; vd=3.33e-4); n=2.17e9u"m^-3", T=(48, 170) .* u"eV"),
    Species(Electron(); n=4.7e9u"m^-3", T=50u"eV");
    B0,
)

# ## Seedless survey
#
# Parallel propagation, `k·d_p ∈ [0.02, 1]`.

dp = 1 / plasma_gyro_ratio((2.53e9 + 2.17e9)u"m^-3", mass(Proton()), B0)
region = (-0.55 - 0.1im, 0.55 + 0.06im)
geom = CartesianSweep(kz=(0.01:0.008:1.0) ./ dp)

sol = solve(DispersionProblem(plasma, region, geom))


# ## Dispersion diagram
using CairoMakie

grow_sol = filter(isgrowing, sol)
fig = Figure(size=(850, 320))
axr = Axis(fig[1, 1]; xlabel="k d_p", ylabel=L"ω_r / Ω_{cp}")
axi = Axis(fig[1, 2]; xlabel="k d_p", ylabel=L"γ / Ω_{cp}")
dispersion_diagram!((axr, axi), grow_sol)
@assert length(grow_sol) == 2 #hide
@assert all(x -> abs(x[1]) < 2e-2, grow_sol) #hide
fig
