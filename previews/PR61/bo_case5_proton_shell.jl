# # Quasi-perpendicular instability — proton shell (BO/ALPS Case 5)
#
# Ion Bernstein / fast-magnetosonic harmonics driven by a spherical proton shell,
# Case 5 of the BO–ALPS comparison (Guo et al. 2026, after Min & Liu 2015):
# a tenuous (10%) 1 keV proton shell at shell speed `v_d = 2 v_A` in a cold
# proton–electron background, propagating at `θ = 89.5°`.
#
# The shell `f_p ∝ exp[−(v − v_d)²/c_p²]` with `v = √(v∥² + v⟂²)` is *not*
# separable in `(v⟂, v∥)`, so it enters as a general analytic `CoupledVDF`.

using VlasovMaxwellDispersion
using DelimitedFiles, Printf
using CairoMakie

# ## Plasma setup
#
# SI parameters from the paper's Table I, normalized to the proton
# gyrofrequency `ωcp`; velocities in units of `c`. The electron mass is the
# table's deliberately heavy `m_e = 10⁻² m_p` — at quasi-perpendicular
# propagation the electron dynamics matter, so the reduced mass must be kept.
#
# The table lists the shell speed as `v_d/c = 0.133`, but the paper's own
# fitted distribution (its Fig. 8i–j) peaks at `|v| = 2 v_A`, matching the
# Min & Liu (2015) setup this case reproduces — we use `v_d = 2 v_A`.

const c0 = 2.99792458e8
qe = 1.602176634e-19; mp = 1.67262192369e-27
eps0 = 8.8541878128e-12; mu0 = 1.25663706212e-6

B0 = 3.28e-8; me = 1.0e-2 * mp
ns = 0.5e5; nc = 4.5e5; ne = 5.0e5           # shell p⁺ / cold p⁺ / e⁻ [m⁻³]

wcp = qe * B0 / mp
vA = B0 / sqrt(mu0 * (ns + nc) * mp)
vd = 2vA / c0                                 # shell speed
cp_ = sqrt(2qe * 1000 / mp) / c0              # shell thermal spread (1 keV)
vthc = sqrt(2qe * 10 / mp) / c0               # cold protons (10 eV)
vthe = sqrt(2qe * 10 / me) / c0               # electrons (10 eV)

Pi2s = ns * qe^2 / (eps0 * mp) / wcp^2
Pi2c = nc * qe^2 / (eps0 * mp) / wcp^2
Pi2e = ne * qe^2 / (eps0 * me) / wcp^2

# `CoupledVDF` takes the unnormalized `f₀(p⟂, p∥)` — the density integral is
# computed internally — over ranges covering the shell (`v_d + 5 c_p`).

hi = vd + 5cp_
shell = CoupledVDF((q, u) -> exp(-(sqrt(q^2 + u^2) - vd)^2 / cp_^2);
    para = (-hi, hi), perp = (0.0, hi))

plasma = (
    NormalizedSpecies(1.0, Pi2s, shell),
    NormalizedSpecies(1.0, Pi2c, Maxwellian(vthc)),
    NormalizedSpecies(-mp / me, Pi2e, Maxwellian(vthe)),
)

# ## Seedless survey
#
# `k` sweeps `k·λ_p ∈ [0.3, 12.5]` (`λ_p = c/ω_pp = v_A/ω_cp`) at `θ = 89.5°`;
# the `ω` box spans the Bernstein harmonic staircase up to `~7.5 ω_cp`.

kunit = c0 / vA                               # k·λ_p → k c/ωcp
region = (0.05 - 0.06im, 7.8 + 0.12im)
geom = AngleSweep(k = (0.3, 12.5) .* kunit, theta = deg2rad(89.5))
sol = solve(GlobalDispersionProblem(plasma, region, geom))

# ## Dispersion diagram

fig = Figure(size = (900, 420))
axr = Axis(fig[1, 1]; xlabel = "k λp", ylabel = "Re ω / ωcp", title = "Proton shell, θ = 89.5°")
axi = Axis(fig[1, 2]; xlabel = "k λp", ylabel = "Im ω / ωcp")
klp(b) = [sqrt(abs2(k)) / kunit for k in b.k]
for b in sol.roots
    x = klp(b)
    p = sortperm(x)
    lines!(axr, x[p], real.(b.omega)[p]; color = (:gray, 0.6), linewidth = 1.5)
    lines!(axi, x[p], imag.(b.omega)[p]; color = (:crimson, 0.8), linewidth = 1.5)
end
ylims!(axr, 0, 8); ylims!(axi, 0, 0.09)
fig
