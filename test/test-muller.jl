@testitem "muller non-finite guards" begin
    using VlasovMaxwellDispersion: muller

    @test isnan(muller(x -> ComplexF64(Inf, Inf), 1.0, 1.1, 1.2))
    @test isnan(muller(x -> ComplexF64(NaN, NaN), 1.0, 1.1, 1.2))

    # flat tail: the quadratic model extrapolates a huge step into overflow;
    # muller should contract overshooting trial steps.
    hits = Ref(0)
    fovf(x) = (v = exp(10x) - 1; isfinite(v) || (hits[] += 1); v)
    r = muller(fovf, -2.0, -1.9, -1.8)
    @test abs(r) < 1.0e-10
    @test hits[] > 0
end
