@testitem "Tensor projection moments" begin
    import VlasovMaxwellDispersion as VMD

    f(v, u) = 1 + 2v + 3u + 4v * u
    vperp = collect(range(0.0, 2.0; length = 5))
    vpara = range(-1.0, 3.0; length = 7)
    F = [f(v, u) for v in vperp, u in vpara]
    fit = fit_grid(BicubicHermite(), vperp, vpara, F)

    integral_power(lo, hi, power) = (hi^(power + 1) - lo^(power + 1)) / (power + 1)
    function reference_moment(m, n)
        iv0, iv1 = integral_power(0.0, 2.0, m), integral_power(0.0, 2.0, m + 1)
        iu0, iu1 = integral_power(-1.0, 3.0, n), integral_power(-1.0, 3.0, n + 1)
        return iv0 * iu0 + 2iv1 * iu0 + 3iv0 * iu1 + 4iv1 * iu1
    end

    for powers in ((1, 0), (3, 0))
        @test VMD.moment(fit, powers) ≈ reference_moment(powers...)
    end
end

@testitem "GridVDF accepts non-tensor projections" begin
    import VlasovMaxwellDispersion as VMD

    struct WrappedMethod end
    struct WrappedProjection{P}
        projection::P
    end
    (p::WrappedProjection)(v, u) = p.projection(v, u)
    VMD.domain(p::WrappedProjection) = VMD.domain(p.projection)
    function VMD.fit_grid(::WrappedMethod, vperp, vpara, F)
        return WrappedProjection(VMD.fit_grid(BicubicHermite(), vperp, vpara, F))
    end

    f(v, u) = exp(-(v^2 + u^2))
    vperp = range(0.0, 4.0; length = 17)
    vpara = collect(range(-4.0, 4.0; length = 25))
    F = [f(v, u) for v in vperp, u in vpara]
    grid = GridVDF(vperp, vpara, F; method = WrappedMethod())

    s = NormalizedSpecies(-1.0, 0.5, grid)
    k = Wavenumber(0.1, 0.4)
    ω = 1.3 + 0.05im
    @test contribution(s, ω, k) ≈
        contribution(NormalizedSpecies(-1.0, 0.5, grid.coupled), ω, k)
end
