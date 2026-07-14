# # Oblique electron ring-beam instability (Guo et al. 2025, Case 3)
#
# Case 3 of the [BO vs ALPS benchmark](https://arxiv.org/abs/2606.14439): an oblique
# instability at `θ = 40°` driven by an electron *ring-beam* — a shifted Maxwellian with a
# parallel drift `v_dz = 0.1c` and a perpendicular ring speed `v_dr = 0.05c` — neutralised by
# a cold-ish Maxwellian electron core. The paper reports that BO and ALPS agree throughout the
# spectrum, so both trustworthy targets collapse onto the single curve reproduced here.

using VlasovMaxwellDispersion
using Printf
using CairoMakie

# ## Plasma setup
#
# SI parameters (paper Table I), normalized to the electron gyrofrequency `|ωce|`; velocities
# in units of `c`. VMD's `GaussianRing` is the literal shifted-Gaussian of Eq. (7),
# `f ∝ exp[-(v∥-v_dz)²/c∥² - (v⊥-v_dr)²/c⊥²]`, with `vth = √(2qT/m)/c = c∥ = c⊥`.
# The total electron density `1e6 m⁻³` is charge-neutralised by an immobile background.

const c0 = 2.99792458e8
qe = 1.602176634e-19; me = 9.1093837015e-31; eps0 = 8.8541878128e-12

B0 = 9.6e-8
n_ring = 1.0e5; n_bg = 9.0e5; n_tot = n_ring + n_bg
T = 51.0                               # eV, isotropic for both electron populations

wce = qe * B0 / me                     # |ωce|, reference gyrofrequency
wpe = sqrt(n_tot * qe^2 / (eps0 * me))
Pi2_ring = n_ring * qe^2 / (eps0 * me) / wce^2
Pi2_bg = n_bg * qe^2 / (eps0 * me) / wce^2
vth = sqrt(2qe * T / me) / c0

vdf_ring = GaussianRing(vth_para = vth, vth_perp = vth, vd = 0.1, vr = 0.05)
plasma = (
    NormalizedSpecies(-1.0, Pi2_ring, vdf_ring),
    NormalizedSpecies(-1.0, Pi2_bg, Maxwellian(vth)),
)

# ## Seedless survey
#
# `k` is swept over `k·λₑ ∈ [0.3, 35]` (`λₑ = c/ωpe`, electron inertial length); `kunit`
# converts to VMD's `k c/ωce` units. Over this range `Λr = k⊥ v_dr/ωce ≲ 4 < 10`, so the
# parabolic-cylinder closure of `GaussianRing` is accurate and `SeparableVDF` is not needed.
# The `ω` box spans `Re ω ∈ [0, 10] ωce` and reaches deep below the real axis so the
# red/green branches stay tracked through their strongly Landau-damped extension to
# `k·λₑ = 35`.

kunit = wpe / wce                      # k·λₑ → k c/ωce
region = (-1.0 - 1.5im, 10.0 + 0.6im)
geom = AngleSweep(k = (0.3, 35.0) .* kunit, theta = deg2rad(40))
sol = solve(GlobalDispersionProblem(plasma, region, geom))

# ## Verification against BO/ALPS
#
# The survey resolves three growing branches. Their peak growth rates and the wavenumbers of
# those peaks are compared with values read from fig. 3(b) of the paper.
kle(b) = [sqrt(abs2(k)) / kunit for k in b.k]   # |k| in units of λₑ⁻¹
peakγ(b) = (g = imag.(b.omega); f = isfinite.(g); any(f) ? maximum(g[f]) : -Inf)

growing = sort(filter(b -> peakγ(b) > 0.05, collect(sol.roots)); by = peakγ, rev = true)
target = [(0.33, 10.0), (0.25, 18.0), (0.13, 10.0)]   # (γ_peak, k·λₑ) from fig. 3(b)
for (b, (γref, kref)) in zip(growing, target)
    x = kle(b); g = imag.(b.omega); j = argmax(replace(g, NaN => -Inf))
    @printf(
        "γ_peak = %.3f (paper ≈ %.2f)  at k·λₑ = %.1f (paper ≈ %.0f)\n",
        g[j], γref, x[j], kref
    )
end

# VMD gives `γ_peak = 0.326, 0.263, 0.133` at `k·λₑ = 10.1, 17.1, 9.0`, matching the paper's
# `≈ 0.33, 0.25, 0.13`. The strongest branch also reproduces the characteristic notch in `γ`
# near `k·λₑ ≈ 4` and the negative-then-positive excursion of the flattening `ω_r → 4 ωce`
# mode — features visible in fig. 3.

# ## Dispersion diagram
#
# Panels mirror fig. 3: (a) real frequency, (b) growth rate, for the three unstable branches.

fig = Figure(size = (760, 340))
axr = Axis(fig[1, 1]; xlabel = "k λₑ", ylabel = "ωr / Ωce", title = "(a)", titlealign = :left)
axi = Axis(fig[1, 2]; xlabel = "k λₑ", ylabel = "γ / Ωce", title = "(b)", titlealign = :left)
for (b, col) in zip(growing, [:red, :royalblue, :green])
    x = kle(b); p = sortperm(x)
    lines!(axr, x[p], real.(b.omega)[p]; color = col, linewidth = 2.5)
    lines!(axi, x[p], imag.(b.omega)[p]; color = col, linewidth = 2.5)
end
hlines!(axi, [0.0]; color = (:black, 0.35), linestyle = :dash)
ylims!(axr, 0, 8); ylims!(axi, -0.4, 0.4)
xlims!(axr, 0, 35); xlims!(axi, 0, 35)
fig
