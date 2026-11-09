# # Ion cyclotron emission — alpha ring-beam (Warwick et al. 2018)
#
# Quasi-perpendicular ion-cyclotron emission (ICE) driven by a fusion-born alpha
# ring-beam. A fast magnetosonic/Bernstein branch propagating at `θ = 89.5°` crosses 
# successive deuteron cyclotron harmonics; each crossing that overlaps the alpha ring 
# resonance goes unstable. The plasma is deuterons + electrons (Maxwellian) plus a 
# hot alpha ring-beam with perpendicular ring speed `v_dr = 0.045c`.
# 
# Reference: Warwick (Fig. 3.8, 2018), Xie (Fig. 3, 2025), Guo (2026, arXiv:2606.14439)

using VlasovMaxwellDispersion
using PlasmaPopulations: gyrofrequency, plasma_frequency, thermal_speed
using VlasovMaxwellDispersion: isgrowing
using Unitful
using CairoMakie

# ## Plasma setup
#
# Deuterons + electrons (Maxwellian) plus the hot alpha ring-beam, a `BiMaxwellian` with
# ring speed `vr_perp`. `ref = Proton()` normalizes to `ωcp = eB₀/mp`.

B0 = 2.1u"T"
α, D, e = Particle(:He, 2), Particle(:D, 1), Electron()
T = 1000u"eV"
vth = thermal_speed(α, T)
plasma = NormalizedPlasma(
    Plasma(
        Population(α, GaussianRing(; vr_perp = 0.045u"c", w_perp = vth, vth_para = vth); n = 1.0e16u"m^-3"),
        Population(D, BiMaxwellian; T, n = 9.98e18u"m^-3"),
        Population(e, BiMaxwellian; T, n = 1.0e19u"m^-3");
        B0,
    );
    ref = Proton(),
)

# ## Growing modes from the seedless survey
#

dp = gyrofrequency(Proton(), B0) / plasma_frequency(Proton(), 1.0e19u"m^-3") # proton inertial length in c/ωcp
kdp = range(3.0, 7.0, length=100) ./ dp
θ = deg2rad(89.5)
region = (1.8 - 0.8im, 4.2 + 0.2im)
geom = AngleSweep(k=kdp, theta=θ)
sol = solve(DispersionProblem(plasma, region, geom))

# The instability lives entirely on the propagating branch, filtered by `γ` above
# a small cutoff that rejects the marginal roots sitting on the flat bands.

grow = [(kdp[i], ω) for b in sol for (i, ω) in enumerate(b.omega) if isfinite(ω) && imag(ω) > 1e-3]
harmonics!(ax) = hlines!(ax, [2.5, 3.0, 3.5, 4.0, 4.5]; color=(:gray, 0.5), linestyle=:dash)

fig = Figure(size=(900, 400))
axr = Axis(fig[1, 1]; xlabel="k d_p", ylabel="ωr / ωcp")
axi = Axis(fig[1, 2]; xlabel="k d_p", ylabel="γ / ωcp")
harmonics!(axr)
scatter!(axr, first.(grow), real.(last.(grow)); color=:crimson)
scatter!(axi, first.(grow), imag.(last.(grow)); color=:crimson)
fig

# ## Clean propagating branch via seeded continuation

kseed = 5.2 / dp
seed = Seed(3.5 + 0.03im, Wavevector(kseed .* sincos(θ)...))
solc = solve(DispersionProblem(plasma, seed, geom))

fig2 = Figure(size=(900, 400))
axr2 = Axis(fig2[1, 1]; xlabel="k d_p", ylabel="ωr / ωcp")
axi2 = Axis(fig2[1, 2]; xlabel="k d_p", ylabel="γ / ωcp")
harmonics!(axr2)
scatter!(axr2, first.(grow), real.(last.(grow)); color=(:crimson, 0.35), markersize=6)  # survey for comparison
lines!(axr2, kdp, real.(solc.omega))
lines!(axi2, kdp, imag.(solc.omega))
fig2
