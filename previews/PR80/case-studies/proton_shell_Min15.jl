# # Quasi-perpendicular instability — proton shell (Min & Liu 2015)
#
# Ion Bernstein / fast-magnetosonic harmonics driven by a spherical proton shell.
# Ref: Min & Liu 2015, also Case 5 in Guo et al. 2026.
# A tenuous (10%) 1 keV proton shell at shell speed `v_d = 2 v_A` in a cold
# proton–electron background, propagating at `θ = 89.5°`.
#
# The shell `f_p ∝ exp[−(v − v_d)²/c_p²]` with `v = √(v∥² + v⟂²)`
# enters as a general analytic function, wrapped in a [`LowRankVDF`](@ref): the
# shell is numerically rank ~12, and factoring it decouples the perpendicular Bessel
# moments (ω-independent, hoisted into the per-`k` plan) from the parallel Landau
# integral, accelerating dispersion-tensor evaluation.

using VlasovMaxwellDispersion
using DelimitedFiles, Printf, Unitful
using CairoMakie

# ## Plasma setup
#
# Normalized to the proton gyrofrequency `ωcp`. The electron mass is the table's
# deliberately heavy `m_e = 10⁻² m_p`.
#
# The shell is a raw `f₀(p⊥, p∥)`, so it is written as a function of the plasma's
# [`scales`](@ref): shell speed `2 v_A`, spread the species' own `vth`.

shell(s) = let vd = 2s.vA, hi = 2s.vA + 5s.vth
    LowRankVDF((q, u) -> exp(-(sqrt(q^2 + u^2) - vd)^2 / s.vth^2); para = (-hi, hi), perp = (0.0, hi))
end

plasma = Plasma(
    Species(Proton(), shell; n = 0.5e5u"m^-3", T = 1000u"eV"),        # 1 keV shell
    Species(Proton(); n = 4.5e5u"m^-3", T = 10u"eV"),                 # cold core
    Species(Particle(z = -1, A = 1.0e-2); n = 5.0e5u"m^-3", T = 10u"eV"),
    B0 = 3.28e-8u"T",
)

# ## Seedless survey
#
# `k` sweeps `k·d_p ∈ [0.3, 12.5]` (`d_p = c/ω_pp = v_A/ω_cp`) at `θ = 89.5°`;
# the `ω` box spans the Bernstein harmonic staircase up to `~7.5 ω_cp`.

dp = scales(plasma).vA                        # d_p = v_A/ω_cp in units of c/ω_cp
region = (0.05 - 0.06im, 7.8 + 0.12im)
geom = AngleSweep(k = range(0.3, 12.5, 128) ./ dp, theta = deg2rad(89.5))
sol = solve(DispersionProblem(plasma, region, geom))

# ## Dispersion diagram

fig = Figure(size = (900, 420))
axr = Axis(fig[1, 1]; xlabel = "k d_p", ylabel = "Re ω / ωcp", title = "Proton shell, θ = 89.5°")
axi = Axis(fig[1, 2]; xlabel = "k d_p", ylabel = "Im ω / ωcp")
kdp(b) = [sqrt(abs2(k)) * dp for k in b.k]
for b in sol.roots
    x = kdp(b)
    lines!(axr, x, real.(b.omega); color = (:gray, 0.6), linewidth = 1.5)
    lines!(axi, x, imag.(b.omega); color = (:crimson, 0.8), linewidth = 1.5)
end
ylims!(axr, 0, 8); ylims!(axi, 0, 0.09)
fig
