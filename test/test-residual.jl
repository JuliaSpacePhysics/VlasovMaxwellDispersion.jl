# Scale-invariant residual
# Fixture: ALPS test_kpar_fast plasma (vA=1e-4) where ‖𝒟‖~1e12, so raw |det 𝒟|
# is ~1e14 at a genuine root (cancellation floor ~‖𝒟‖³ε) — every threshold below
# only works because sol.resid is Hadamard-normalized (|det|/∏‖rowᵢ‖ ~ ε at a root).

@testitem "local solve residual is scale-invariant" begin
    import Roots: Order1
    vA = 1.0e-4
    me = 5.44662e-4
    plasma = (
        NormalizedSpecies(1.0, 1 / vA^2, Maxwellian(vA)),
        NormalizedSpecies(-1 / me, 1 / (me * vA^2), Maxwellian(vA / sqrt(me))),
    )
    k = Wavenumber(0.01 / vA, 0.01 / vA)
    prob = LocalDispersionProblem(plasma, k, 9.9881e-3 - 2.3132e-7im)  # ALPS fast-wave root

    for alg in (Muller(), Order1())
        sol = solve(prob, alg)
        @test sol.retcode == :Success
        @test sol.resid < 1.0e-10
    end
end

@testitem "global and branch solves report scale-invariant residuals" begin
    vA = 1.0e-4
    me = 5.44662e-4
    plasma = (
        NormalizedSpecies(1.0, 1 / vA^2, Maxwellian(vA)),
        NormalizedSpecies(-1 / me, 1 / (me * vA^2), Maxwellian(vA / sqrt(me))),
    )
    k = Wavenumber(0.01 / vA, 0.01 / vA)
    ωref = 9.9881e-3 - 2.3132e-7im

    gsol = solve(GlobalDispersionProblem(plasma, k, (0.008 - 0.001im, 0.012 + 0.001im)), GRPF(; tol=1.0e-4))
    @test gsol.retcode == :Success
    @test length(gsol.resid) == length(gsol.omega)
    # GRPF roots are mesh-accurate (~tol), not polished, hence the loose bound.
    @test all(0 .<= gsol.resid .< 1.0e-2)

    bsol = solve(BranchProblem(plasma, [k], ωref))
    @test bsol.retcode == :Success
    @test only(bsol.resid) < 1.0e-10
end
