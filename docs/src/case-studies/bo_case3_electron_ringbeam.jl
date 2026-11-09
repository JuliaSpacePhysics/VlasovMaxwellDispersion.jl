# # Oblique electron ring-beam instability
#
# Oblique instability at `θ = 40°` driven by an electron *ring-beam* — 
# a shifted Maxwellian with a parallel drift `v_dz = 0.1c` and a perpendicular 
# ring speed `v_dr = 0.05c` — neutralised by a cold-ish Maxwellian electron core.
# 
# Reference: Guo (Case 3, arXiv:2606.14439)

using VlasovMaxwellDispersion
using PlasmaPopulations: gyrofrequency, plasma_frequency, thermal_speed
using VlasovMaxwellDispersion: isgrowing
using Unitful, Printf
using CairoMakie

# ## Plasma setup
#
# The ring beam is a `GaussianRing` with ring speed `vr_perp`, the literal shifted Gaussian
# `f ∝ exp[-(v∥-vd_para)²/vth² - (v⊥-vr_perp)²/w_perp²]`, with `w_perp = vth = √(2qT/m)`.

B0 = 9.6e-8u"T"
e = Electron()
vth = thermal_speed(e, 51u"eV")
plasma = Plasma(
    Population(e, GaussianRing(; vr_perp = 0.05u"c", w_perp = vth, vth_para = vth, vd_para = 0.1u"c"); n = 1.0e5u"m^-3"),
    Population(e, BiMaxwellian(; vth); n = 9.0e5u"m^-3");
    B0,
)

# ## Seedless survey

de = abs(gyrofrequency(e, B0)) / plasma_frequency(e, 1.0e6u"m^-3")  # c/ωpe in units of c/ωce
ks = range(0.5, 35.0, 128) ./ de
region = (-1.0 - 1.5im, 10.0 + 0.6im)
geom = AngleSweep(k=ks, theta=deg2rad(40))
sol = link(solve(DispersionProblem(NormalizedPlasma(plasma), region, geom), AAA()))

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