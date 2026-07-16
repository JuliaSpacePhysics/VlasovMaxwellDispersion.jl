@testitem "Track follows vacuum light branch" begin
    using VlasovMaxwellDispersion: ReturnCode
    plasma = NormalizedSpecies(0.0, 0.0, ColdVDF())
    kzs = range(0.5, 1.0; length=12)
    ks = Wavenumber.(0.0, kzs)
    sol = solve(DispersionProblem(plasma, 0.5 + 0im, ks)) # default to Continuation with order=3
    @test sol.omega ≈ kzs rtol = 1.0e-5
    @testset "track predictor order: degree 1 (secant) through 3 all hold the branch" begin
        for order in 1:3
            sol = solve(DispersionProblem(plasma, 0.5 + 0im, ks), Continuation(; order))
            @test sol.retcode == ReturnCode.Success
            @test sol.omega ≈ kzs rtol = 1.0e-6
        end
    end

    ks = (Wavenumber(0.0, kz) for kz in (0.5, 0.6, 0.7))
    roots = solve(DispersionProblem(plasma, 0.5 + 0im, ks)).omega
    @test all(isfinite, roots)

    @testset "track accepts a swept geometry as the k path" begin
        sol = solve(DispersionProblem(plasma, 0.5 + 0im, CartesianSweep(kz=(0.5, 1.0))))
        @test sol.retcode == ReturnCode.Success
        @test sol.omega ≈ range(0.5, 1.0; length=61) atol = 1.0e-6
        @test_broken solve(
            DispersionProblem(plasma, 0.5 + 0im, CartesianSweep(kx=(0.0, 0.1), kz=(0.5, 1.0)))
        )
    end
end


@testitem "track subdivides across Gary84 v0=10 branch transition" begin
    mp_me = 1836.15267343
    vA_c = 1.0e-4
    nm, nb = 0.99, 0.01
    vm_c = vA_c * sqrt(1 / (2nm))
    vth_m = sqrt(2) * vm_c
    vth_b = sqrt(10) * vth_m
    vth_e = vth_m * sqrt(mp_me)
    pi2_i = 1 / vA_c^2
    v0 = 10 * vm_c
    v0m = -nb * v0 / (nm + nb)
    v0b = v0m + v0

    plasma = (
        NormalizedSpecies(1.0, nm * pi2_i, Maxwellian(; vth_para=vth_m, vth_perp=vth_m, vd=v0m)),
        NormalizedSpecies(1.0, nb * pi2_i, Maxwellian(; vth_para=vth_b, vth_perp=vth_b, vd=v0b)),
        NormalizedSpecies(-mp_me, pi2_i * mp_me, Maxwellian(; vth_para=vth_e, vth_perp=vth_e)),
    )
    ks = [Wavenumber(0.0, ka / vm_c) for ka in 0.1:-0.001:0.01]

    # This grid is too coarse through the transition to keep the branch: with
    # subdivision off, muller contracts its overflowing trial steps and wanders
    # to a distant root (the electron plasma branch) while still "converging".
    baseline = solve(DispersionProblem(plasma, 0.08im, ks), Continuation(; maxsubdiv=0))
    @test !isapprox(baseline.omega[end], 0.009597 + 0im; rtol=5.0e-3, atol=5.0e-5)

    roots = solve(DispersionProblem(plasma, 0.08im, ks), Continuation()).omega

    @test all(isfinite, roots)
    @test roots[58] ≈ 0.015838 + 0.000212im rtol=5.0e-3 atol=5.0e-5 # ka=0.043
    @test roots[end] ≈ 0.009597 + 0im rtol=5.0e-3 atol=5.0e-5       # ka=0.010
end
