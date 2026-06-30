@testitem "arc-length track follows vacuum light branch" begin
    using VlasovMaxwellDispersion

    plasma = NormalizedSpecies(0.0, 0.0, ColdVDF())
    ks = [Wavenumber(0.0, kz) for kz in range(0.5, 1.0; length=6)]
    roots = solve(BranchProblem(plasma, ks, 0.5 + 0im), ArcLength()).omega

    @test all(isfinite, roots)
    @test maximum(abs.(roots .- [k.kz for k in ks])) < 1e-5
end

@testitem "arc-length track accepts unsized iterables" begin
    plasma = NormalizedSpecies(0.0, 0.0, ColdVDF())
    ks = (Wavenumber(0.0, kz) for kz in (0.5, 0.6, 0.7))
    roots = solve(BranchProblem(plasma, ks, 0.5 + 0im), ArcLength()).omega

    @test length(roots) == 3
    @test all(isfinite, roots)
end


@testitem "track fallback crosses Gary84 v0=10 branch transition" begin
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

    # Without fallback the tracker loses the branch at the transition — historically as NaN;
    # since muller contracts overflowing trial steps it instead wanders to a distant root
    # (the electron plasma branch). Either failure mode is what the jump fallback must catch.
    baseline = solve(BranchProblem(plasma, ks, 0.08im),
                     ArcLength(; fallback=false, atol=1.0e-10, maxiter=300)).omega
    @test !isapprox(baseline[end], 0.009597 + 0im; rtol=5.0e-3, atol=5.0e-5)

    roots = solve(BranchProblem(plasma, ks, 0.08im), ArcLength(; atol=1.0e-10, maxiter=300)).omega

    @test all(isfinite, roots)
    @test roots[58] ≈ 0.015838 + 0.000212im rtol=5.0e-3 atol=5.0e-5 # ka=0.043
    @test roots[end] ≈ 0.009597 + 0im rtol=5.0e-3 atol=5.0e-5       # ka=0.010
end