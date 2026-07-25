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
using CairoMakie

# ## Plasma setup
#
# Normalized to the proton gyrofrequency `ωcp`; drift speeds are in `c`. Current closure
# is satisfied by the two proton drifts (`n_c v_c + n_b v_b ≈ 0`); electrons are isotropic
# and at rest.

plasma = Plasma(
    Species(Proton(), Maxwellian(vd=-2.87e-4); n=2.53e9u"m^-3", T=(20, 100) .* u"eV"),
    Species(Proton(), Maxwellian(vd=3.33e-4); n=2.17e9u"m^-3", T=(48, 170) .* u"eV"),
    Species(Particle(z=-1, A=5.447e-4); n=4.7e9u"m^-3", T=50u"eV"),
    B0=7.5e-7u"T",
)

# ## Seedless survey
#
# Parallel propagation, `k·λₚ ∈ [0.02, 1]`. `scales(plasma).di` is `λₚ = c/ωpp` over
# *both* proton populations — the paper's total-density convention.

λp = scales(plasma).di
region = (-0.55 - 0.1im, 0.55 + 0.06im)
geom = CartesianSweep(kz=(0.01:0.008:1.0) ./ λp)

sol = solve(DispersionProblem(plasma, region, geom))


# ## Dispersion diagram
grow_sol = filter(isgrowing, sol)
fig = Figure(size=(850, 320))
axr = Axis(fig[1, 1]; xlabel="k λₚ", ylabel=L"ω_r / Ω_{cp}")
axi = Axis(fig[1, 2]; xlabel="k λₚ", ylabel=L"γ / Ω_{cp}")
dispersion_diagram!((axr, axi), grow_sol)
@assert length(grow_sol) == 2 #hide
@assert all(x -> abs(x[1]) < 2e-2, grow_sol) #hide
fig
