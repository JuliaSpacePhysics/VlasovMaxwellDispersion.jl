# # Ion cyclotron waves from a proton core-beam
# 
# Counter-propagating ion cyclotron waves excited by a proton core-beam
# distribution. Both proton populations are drifting bi-Maxwellians with
# `T⊥ > T∥`; the combined anisotropy and beam drift destabilize forward- and
# backward-propagating ICWs with unequal growth rates.

# Reference: Case 6 in Guo et al. (2026, arXiv:2606.14439)

using VlasovMaxwellDispersion
using PlasmaPopulations: gyrofrequency, plasma_frequency
using VlasovMaxwellDispersion: isgrowing
using Unitful

# ## Plasma setup
#
B0 = 7.5e-7u"T"
p, e = Proton(), Electron()
plasma = Plasma(
    Population(p, BiMaxwellian; T_perp = 100u"eV", T_para = 20u"eV", vd_para = -2.87e-4u"c", n = 2.53e9u"m^-3"),
    Population(p, BiMaxwellian; T_perp = 170u"eV", T_para = 48u"eV", vd_para = 3.33e-4u"c", n = 2.17e9u"m^-3"),
    Population(e, BiMaxwellian; T = 50u"eV", n = 4.7e9u"m^-3");
    B0,
)

# ## Seedless survey
#
# Parallel propagation, `k·d_p ∈ [0.02, 1]`.

dp = gyrofrequency(p, B0) / plasma_frequency(p, (2.53e9 + 2.17e9)u"m^-3")  # c/ωpp in units of c/ωcp
region = (-0.55 - 0.1im, 0.55 + 0.06im)
geom = CartesianSweep(k_para=(0.01:0.008:1.0) ./ dp)

sol = link(solve(DispersionProblem(NormalizedPlasma(plasma), region, geom), AAA()))


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
