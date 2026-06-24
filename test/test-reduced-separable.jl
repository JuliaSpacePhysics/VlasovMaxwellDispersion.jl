@testitem "generic analytic hilbert reproduces Z" begin
    using VlasovMaxwellDispersion: hilbert, Z
    g(v) = exp(-v^2) / sqrt(pi)                      # ∫ g/(v-ζ) = Z(ζ), Im ζ>0
    for ζ in (1.5 + 0.8im, 0.3 + 0.05im, 1.2 - 0.3im) # upper, near-real, lower (Landau)
        @test abs(hilbert(g, ζ; lower = -20.0, upper = 20.0) - Z(ζ)) < 1.0e-10
    end
end

@testitem "reduced SeparableVDF(Gaussian) electrostatic root matches Maxwellian" begin
    using VlasovMaxwellDispersion
    vth, Pi2, kz = 1.0, 1.0, 0.7
    k = Wavenumber(0.0, kz)
    f = SeparableVDF(u -> exp(-(u / vth)^2) / (sqrt(pi) * vth); lower = -12.0, upper = 12.0)
    ωf = solve(LocalDispersionProblem(NormalizedSpecies(-1.0, Pi2, f), k, 1.2 - 0.1im)).omega
    ωm = solve(LocalDispersionProblem(NormalizedSpecies(-1.0, Pi2, Maxwellian(vth)), k, 1.2 - 0.1im)).omega
    @test abs(ωf - ωm) < 1.0e-8
end

@testitem "reduced SeparableVDF bump-on-tail is unstable (Im ω > 0)" begin
    using VlasovMaxwellDispersion
    vth, vb, wb, nb = 1.0, 4.0, 0.5, 0.06
    fbump(u) = (1 - nb) * exp(-(u / vth)^2) / (sqrt(pi) * vth) +
        nb * exp(-((u - vb) / wb)^2) / (sqrt(pi) * wb)
    p = NormalizedSpecies(-1.0, 1.0, SeparableVDF(fbump; lower = -12.0, upper = 14.0))
    for kz in (0.2, 0.25, 0.3)                        # phase speed on the bump front
        ω = solve(LocalDispersionProblem(p, Wavenumber(0.0, kz), kz * vb + 0.05im)).omega
        @test imag(ω) > 0
    end
end

@testitem "reduced SeparableVDF guards against non-electrostatic use" begin
    using VlasovMaxwellDispersion
    using VlasovMaxwellDispersion: contribution
    f = SeparableVDF(u -> exp(-u^2) / sqrt(pi); lower = -10.0, upper = 10.0)
    s = NormalizedSpecies(-1.0, 1.0, f)
    @test_throws ArgumentError contribution(s, 1.0 + 0im, Wavenumber(0.3, 0.5))  # kperp≠0
end
