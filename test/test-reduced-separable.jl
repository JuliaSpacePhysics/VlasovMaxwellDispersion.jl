@testitem "generic analytic Landau integral reproduces Z" begin
    using VlasovMaxwellDispersion: landau, Z
    g(v) = exp(-v^2) / sqrt(pi)                      # ∫ g/(v-ζ) = Z(ζ), Im ζ>0
    L, U = -30.0, 30.0
    # upper, near-real, real, lower
    for ζ in (1.5 + 0.8im, 0.3 + 0.05im, 0.3, 0.0, 1.2 - 0.3im, 2.0 - 5.0im)
        @test landau(g, ζ, L, U) ≈ Z(ζ)
        @test landau(g, [ζ], L, U) ≈ [Z(ζ)]    # vector path, single pole
    end
    # Far and near-but-unpeelable causal poles avoid `Inf * 0`.
    # Damped pole still carries its Landau residue.
    ζs = [80.0im, 30.0im, 1.2 - 0.3im]
    @test !isfinite(g(30.0im))
    @test landau(g, ζs, L, U) ≈ Z.(ζs)
    @test isnan(Z(1.0 - 30.0im))
    @test isnan(landau(g, 1.0 - 30.0im, L, U))
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
