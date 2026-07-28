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
# `GaussianRing` is the literal shifted-Gaussian, `f ∝ exp[-(v∥-v_dz)²/c∥² - (v⊥-v_dr)²/c⊥²]`,
# with `c∥ = c⊥ = vth = √(2qT/m)` from `T` and the drifts in `c`.

B0 = 9.6e-8u"T"
plasma = Plasma(
    Species(Electron(), sc -> GaussianRing(sc; vd=0.1, vr=0.05); n=1.0e5u"m^-3", T=51u"eV"),
    Species(Electron(); n=9.0e5u"m^-3", T=51u"eV");
    B0,
)

# ## Seedless survey

de = 1 / plasma_gyro_ratio(1.0e6u"m^-3", mass(Electron()), B0)
ks = range(0.5, 35.0, 128) ./ de
region = (-1.0 - 1.5im, 10.0 + 0.6im)
geom = AngleSweep(k=ks, theta=deg2rad(40))
sol = solve(DispersionProblem(plasma, region, geom))

# The survey resolves three growing branches. Here we compare their peak growth rates and the wavenumbers of those peaks.

kde(b) = [sqrt(abs2(k)) * de for k in b.k]   # |k| in units of d_e⁻¹

growing = filter(x -> isgrowing(x, 0.05), sol)
for b in growing
    x = kde(b);
    g = imag.(b.omega);
    j = argmax(replace(g, NaN => -Inf))
    @printf("γ_peak = %.3f  at k·d_e = %.1f\n", g[j], x[j])
end

# ## Dispersion diagram
dispersion_diagram(growing)