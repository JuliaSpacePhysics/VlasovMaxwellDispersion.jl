# det 𝒟 is complex-valued: SimpleNonlinearSolve's Broyden/DFSane compare residual magnitudes
# via `isless(::Complex,…)` and break; SimpleHalley (FiniteDiff AD) works.
@testitem "SciML interop" begin
    import SimpleNonlinearSolve as SNS

    p = NormalizedPopulation(1.0, 1.0, BiMaxwellian(; vth = 1.0))
    prob = DispersionProblem(p, 0.6, Wavevector(0.01, 0.5))
    sd = solve(prob, SNS.SimpleHalley())
    @test sd.retcode == ReturnCode.Success
    @test sd.omega ≈ solve(prob, Muller()).omega rtol = 1e-6
end
