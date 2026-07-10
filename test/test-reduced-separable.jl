@testitem "generic analytic hilbert reproduces Z" begin
    using VlasovMaxwellDispersion: plan_landau, Z
    g(v) = exp(-v^2) / sqrt(pi)                      # ∫ g/(v-ζ) = Z(ζ), Im ζ>0
    # upper, near-real, lower (Landau)
    for ζ in (1.5 + 0.8im, 0.3 + 0.05im, 1.2 - 0.3im, 2.0 - 5.0im, 3.0 - 9.0im)
        @test plan_landau((-30.0, 30.0), ζ)(g) ≈ Z(ζ) rtol = 1.0e-10
    end
end

@testitem "ReducedVDF(Gaussian) electrostatic root matches Maxwellian" begin
    using VlasovMaxwellDispersion
    vth, Pi2, kz = 1.0, 1.0, 0.7
    k = Wavenumber(0.0, kz)
    f = ReducedVDF(u -> exp(-(u / vth)^2) / (sqrt(pi) * vth); para = (-12.0, 12.0))
    ωf = solve(DispersionProblem(NormalizedSpecies(-1.0, Pi2, f), 1.2 - 0.1im, k)).omega
    ωm = solve(DispersionProblem(NormalizedSpecies(-1.0, Pi2, Maxwellian(vth)), 1.2 - 0.1im, k)).omega
    @test abs(ωf - ωm) < 1.0e-8
end

@testitem "ReducedVDF bump-on-tail is unstable (Im ω > 0)" begin
    using VlasovMaxwellDispersion
    vth, vb, wb, nb = 1.0, 4.0, 0.5, 0.06
    fbump(u) = (1 - nb) * exp(-(u / vth)^2) / (sqrt(pi) * vth) +
        nb * exp(-((u - vb) / wb)^2) / (sqrt(pi) * wb)
    p = NormalizedSpecies(-1.0, 1.0, ReducedVDF(fbump; para = (-12.0, 14.0)))
    for kz in (0.2, 0.25, 0.3)                        # phase speed on the bump front
        ω = solve(DispersionProblem(p, kz * vb + 0.05im, Wavenumber(0.0, kz))).omega
        @test imag(ω) > 0
    end
end

@testitem "ReducedVDF guards against non-electrostatic use" begin
    using VlasovMaxwellDispersion
    using VlasovMaxwellDispersion: contribution
    f = ReducedVDF(u -> exp(-u^2) / sqrt(pi); para = (-10.0, 10.0))
    s = NormalizedSpecies(-1.0, 1.0, f)
    @test_throws ArgumentError contribution(s, 1.0 + 0im, Wavenumber(0.3, 0.5))  # kperp≠0
end
