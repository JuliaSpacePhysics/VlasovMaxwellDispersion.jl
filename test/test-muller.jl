@testitem "muller non-finite guards" begin
    using VlasovMaxwellDispersion: muller

    @test isnan(muller(x -> ComplexF64(Inf, Inf), 1.0, 1.1, 1.2))
    @test isnan(muller(x -> ComplexF64(NaN, NaN), 1.0, 1.1, 1.2))

    # flat tail: quadratic model extrapolates huge steps into overflow;
    # muller should contract overshooting trial steps.
    hits = Ref(0)
    fovf(x) = (v = exp(10x) - 1; isfinite(v) || (hits[] += 1); v)
    r = muller(fovf, -2.0, -1.9, -1.8)
    @test abs(r) < 1.0e-10
    @test hits[] > 0
end

@testitem "Polishes |ω|≪1 roots (relative seed spread)" begin
    # Dense cold e-p plasma: EMIC/Alfvén root near ω≈8.6e-5 (≪ 1).
    mp_me = 1836.15
    plasma = (
        NormalizedSpecies(-1.0, 100.0, ColdVDF()),
        NormalizedSpecies(1 / mp_me, 100.0 / mp_me, ColdVDF()),
    )
    k = Wavenumber(0.0, 0.04)
    sol = solve(DispersionProblem(plasma, 8.5e-5 + 0.0im, k))
    @test sol.retcode == :Success
    @test abs(sol.omega - 8.5e-5) < 1.0e-5   # stayed on low-ω branch
    @test sol.resid < 1.0e-8                 # polished to genuine root
end
