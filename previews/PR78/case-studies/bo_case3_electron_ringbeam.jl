# # Oblique electron ring-beam instability
#
# Oblique instability at `θ = 40°` driven by an electron *ring-beam* — 
# a shifted Maxwellian with a parallel drift `v_dz = 0.1c` and a perpendicular 
# ring speed `v_dr = 0.05c` — neutralised by a cold-ish Maxwellian electron core.
# 
# Reference: Guo (Case 3, arXiv:2606.14439)

using VlasovMaxwellDispersion
using VlasovMaxwellDispersion: isgrowing
using Unitful, Printf
using CairoMakie

# ## Plasma setup
#
# Normalized to the electron gyrofrequency `|ωce|`. VMD's `GaussianRing` is the literal
# shifted-Gaussian of Eq. (7), `f ∝ exp[-(v∥-v_dz)²/c∥² - (v⊥-v_dr)²/c⊥²]`, with
# `c∥ = c⊥ = vth = √(2qT/m)` from `T` and the drifts in `c`. The total electron density
# `1e6 m⁻³` is charge-neutralised by an immobile background.

plasma = Plasma(
    Species(Electron(), GaussianRing(vd=0.1, vr=0.05); n=1.0e5u"m^-3", T=51u"eV"),
    Species(Electron(); n=9.0e5u"m^-3", T=51u"eV"),
    B0=9.6e-8u"T",
)

# ## Seedless survey
#
# `k` is swept over `k·λₑ ∈ [0.3, 35]`; `scales(plasma).de` is the electron inertial
# length `λₑ = c/ωpe` (at the *total* electron density) in VMD's `c/ωce` units.
# Over this range `Λr = k⊥ v_dr/ωce ≲ 4 < 10`, so the parabolic-cylinder closure of
# `GaussianRing` is accurate and `SeparableVDF` is not needed.
# The `ω` box spans `Re ω ∈ [0, 10] ωce` and reaches deep below the real axis so the
# red/green branches stay tracked through their strongly Landau-damped extension to
# `k·λₑ = 35`.

λe = scales(plasma).de
region = (-1.0 - 1.5im, 10.0 + 0.6im)
geom = AngleSweep(k=range(0.5, 35.0, 128) ./ λe, theta=40u"°")
sol = solve(DispersionProblem(plasma, region, geom))

# The survey resolves three growing branches. Here we compare their peak growth rates and the wavenumbers of those peaks.

kle(b) = [sqrt(abs2(k)) * λe for k in b.k]   # |k| in units of λₑ⁻¹

growing = filter(x -> isgrowing(x, 0.05), sol)
for b in growing
    x = kle(b);
    g = imag.(b.omega);
    j = argmax(replace(g, NaN => -Inf))
    @printf("γ_peak = %.3f  at k·λₑ = %.1f\n", g[j], x[j])
end

# ## Dispersion diagram
dispersion_diagram(growing)