@testitem "holomorphic AD rules for special functions vs finite differences" begin
    using VlasovMaxwellDispersion: _dwrt, _val_dwrt
    using SpecialFunctions: erf, erfc, erfcx, erfi, dawson
    using Gamma: gamma
    z = 0.7 - 0.3im
    h = 1.0e-6
    for f in (erf, erfc, erfcx, erfi, dawson, gamma)
        fd = (f(z + h) - f(z - h)) / 2h
        v, d = _val_dwrt(f, z)
        @test v == f(z)
        @test d ≈ fd rtol = 1.0e-8
        @test _dwrt(f, z) == d
    end
    @test _dwrt(erfc, 0.5) isa Float64    # real argument keeps the real fast path
end

@testitem "holomorphic AD composes through erfc-built expressions" begin
    using VlasovMaxwellDispersion: _dwrt, _grad2
    using SpecialFunctions: erfc, dawson
    using Gamma: gamma
    h = 1.0e-6
    g(u) = exp(-u^2) * erfc(-u)
    u = 1.1 - 0.4im
    @test _dwrt(g, u) ≈ (g(u + h) - g(u - h)) / 2h rtol = 1.0e-8
    f2(x, y) = erfc(x * y) * gamma(y) + dawson(x)
    x, y = 0.4 - 0.2im, 1.3 + 0.1im
    dx, dy = _grad2(f2, x, y)
    @test dx ≈ (f2(x + h, y) - f2(x - h, y)) / 2h rtol = 1.0e-8
    @test dy ≈ (f2(x, y + h) - f2(x, y - h)) / 2h rtol = 1.0e-8
end

@testitem "ReducedVDF erfc-skewed f₀: auto-df root matches hand-supplied df" begin
    using VlasovMaxwellDispersion
    using SpecialFunctions: erfc
    fsk(u) = exp(-u^2) * erfc(-u) / sqrt(π)                                        # ∫ = 1 by u → -u symmetry
    dfsk(u) = (-2u * erfc(-u) + 2 * exp(-u^2) / sqrt(π)) * exp(-u^2) / sqrt(π)
    auto = ReducedVDF(fsk; para = (-10.0, 12.0))
    manual = ReducedVDF(fsk; para = (-10.0, 12.0), df = dfsk)
    k = Wavenumber(0.0, 0.5)
    ωa = solve(DispersionProblem(NormalizedSpecies(-1.0, 1.0, auto), 1.3 - 0.1im, k)).omega
    ωm = solve(DispersionProblem(NormalizedSpecies(-1.0, 1.0, manual), 1.3 - 0.1im, k)).omega
    @test abs(ωa - ωm) < 1.0e-12
end
