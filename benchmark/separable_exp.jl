# Experiment harness for the SeparableVDF fused-quadrature + Bessel-ladder work.
# Run: julia --project=benchmark benchmark/separable_exp.jl
using VlasovMaxwellDispersion
using VlasovMaxwellDispersion: contribution
using Chairmarks

vthp, vthq = 0.9, 1.2
mx = Maxwellian(vth_para = vthp, vth_perp = vthq)
sep = prepare(SeparableVDF(mx; para = (-14vthp, 14vthp), perp = 14vthq))

# Non-Gaussian parallel (kappa-like) × Gaussian perp — no closed form, exercises quadrature.
fpar(u) = (1 + u^2 / 3)^(-2)
sepk = prepare(SeparableVDF(v -> exp(-v^2) / pi, fpar; para = (-30.0, 30.0), perp = 10.0))

cases = [
    ("sep/kp=0.1", NormalizedSpecies(-1.0, 0.5, sep), 1.3 - 0.05im, Wavenumber(0.1, 0.4)),
    ("sep/kp=1.0", NormalizedSpecies(-1.0, 0.5, sep), 1.3 - 0.05im, Wavenumber(1.0, 0.4)),
    ("sep/kp=2.5", NormalizedSpecies(-1.0, 0.5, sep), 1.3 - 0.05im, Wavenumber(2.5, 0.4)),
    ("kappa/kp=1.0", NormalizedSpecies(-1.0, 1.0, sepk), 1.2 - 0.05im, Wavenumber(1.0, 0.4)),
    ("kappa/kp=3.0", NormalizedSpecies(-1.0, 1.0, sepk), 1.2 - 0.05im, Wavenumber(3.0, 0.4)),
]

# Correctness vs Maxwellian closed form (Gaussian cases only).
println("# correctness (vs Maxwellian closed form)")
for (Ω, Pi2, ω, kz, kp) in (
        (-1.0, 0.5, 1.3 - 0.05im, 0.4, 0.3),
        (-1.0, 0.5, 0.7 + 0.02im, 0.25, 0.6),
        (2.0, 0.8, 1.1 - 0.1im, 0.5, 0.2),
        (-1.0, 0.5, 1.3 - 0.05im, 0.4, 2.5),
    )
    k = Wavenumber(kp, kz)
    χs = contribution(NormalizedSpecies(Ω, Pi2, sep), ω, k)
    χm = contribution(NormalizedSpecies(Ω, Pi2, mx), ω, k)
    rel = maximum(abs, χs .- χm) / maximum(abs, χm)
    println("  Ω=$Ω kp=$kp  relerr=$rel")
end

println("\n# timings")
for (name, s, ω, k) in cases
    b = @b contribution($s, $ω, $k)
    println("  $name: ", b)
end
