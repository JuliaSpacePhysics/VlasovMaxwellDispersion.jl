# # Oblique electron ring-beam instability
#
# Oblique instability at `θ = 40°` driven by an electron *ring-beam* — 
# a shifted Maxwellian with a parallel drift `v_dz = 0.1c` and a perpendicular 
# ring speed `v_dr = 0.05c` — neutralised by a cold-ish Maxwellian electron core.
# 
# Reference: Guo (Case 3, arXiv:2606.14439)

using VlasovMaxwellDispersion
using VlasovMaxwellDispersion: isgrowing
using Printf
using CairoMakie

# ## Plasma setup
#
# SI parameters, normalized to the electron gyrofrequency `|ωce|`; velocities
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

vdf_ring = GaussianRing(vth_para=vth, vth_perp=vth, vd=0.1, vr=0.05)
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
geom = AngleSweep(k=range(0.3, 35.0, 100) .* kunit, theta=deg2rad(40))
sol = solve(DispersionProblem(plasma, region, geom))

# The survey resolves three growing branches. Here we compare their peak growth rates and the wavenumbers of those peaks.

kle(b) = [sqrt(abs2(k)) / kunit for k in b.k]   # |k| in units of λₑ⁻¹

growing = filter(x -> isgrowing(x, 0.05), sol)
for b in growing
    x = kle(b);
    g = imag.(b.omega);
    j = argmax(replace(g, NaN => -Inf))
    @printf("γ_peak = %.3f  at k·λₑ = %.1f\n", g[j], x[j])
end

# ## Dispersion diagram
dispersion_diagram(growing; figure=(size=(300, 500),))