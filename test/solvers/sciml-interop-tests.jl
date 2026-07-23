# det 𝒟 is holomorphic but complex-valued, so most NonlinearSolve solvers break: 
# the bracketing methods need a real sign change except `BracketingNonlinearSolve.Muller`
# and common SimpleNonlinearSolve method Broyden/DFSane compare residual magnitudes via `isless(::Complex,…)`
# except SimpleHalley (FiniteDiff AD).

@testitem "SciML interop" begin
    import SimpleNonlinearSolve as SNS
    import BracketingNonlinearSolve as BNS

    p = NormalizedSpecies(Omega=1.0, Pi2=1.0, vdf=Maxwellian(1.0))
    k = Wavenumber(kperp=0.01, kz=0.5)
    prob = DispersionProblem(p, 0.6, k)
    f = prob.f
    ref = solve(prob).omega
    @test abs(f(ref)) < 1e-8

    h = 1e-3
    ip = SNS.IntervalNonlinearProblem((ω, _) -> f(ω), (prob.target[] - h, prob.target[] + h))
    sm = solve(ip, BNS.Muller(prob.target[] + h * im))
    @test sm.retcode == SNS.ReturnCode.Success
    @test sm.u ≈ ref rtol=1e-6

    # Initial-guess path
    np = SNS.NonlinearProblem((ω, _) -> f(ω), complex(prob.target[]))
    sh = solve(np, SNS.SimpleHalley())
    @test sh.u ≈ ref rtol=1e-6

    # SciMLBaseExt with SciML alg
    sd = solve(prob, SNS.SimpleHalley())
    @test sd.retcode == ReturnCode.Success
    @test sd.omega ≈ ref rtol=1e-6
    @test sd.stats.nevals > 0
end
