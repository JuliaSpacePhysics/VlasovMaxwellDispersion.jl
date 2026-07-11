@testitem "GRPF deflates the ω=0 light-term pole" begin
    using VlasovMaxwellDispersion: NormalizedSpecies, Wavenumber, ColdVDF,
        GlobalDispersionProblem, GRPF, solve, dispersion_tensor,
        wave_dispersion_tensor
    using RootsAndPoles
    using LinearAlgebra

    # Cold e-p plasma dense enough (low vA) that its parallel modes sit near ω=0.
    mp_me = 1836.15
    plasma = (
        NormalizedSpecies(-1.0, 100.0, ColdVDF()),
        NormalizedSpecies(1 / mp_me, 100.0 / mp_me, ColdVDF()),
    )
    k = Wavenumber(0.0, 2.0)
    region = (-0.08 - 0.03im, 0.08 + 0.03im)   # straddles ω=0
    tol = 0.004

    # det𝒟 blows up as ω→0 (light-term pole); deflation clears it so det(ω̃²𝒟) stays finite.
    detD(w) = det(dispersion_tensor(plasma, w, k))
    detD_wave(w) = det(wave_dispersion_tensor(plasma, w, k))
    @test abs(detD(1.0e-3)) > 1.0e3 * abs(detD(1.0e-2))
    @test isnan(detD(0))
    @test isfinite(detD_wave(1.0e-10))
    @test detD_wave(0) ≈ detD_wave(1.0e-10)

    # Bug: raw det𝒟's ω=0 pole winding cancels nearby roots
    # GRPF reports a spurious pole near origin
    oc = rectangulardomain(ComplexF64(region[1]), ComplexF64(region[2]), tol)
    _, pD = grpf(detD, oc, GRPFParams(9000, tol, false))
    @test any(p -> abs(p) < 0.02, pD)

    # Fix: the deflated global solver is pole-free at the origin
    sol = solve(GlobalDispersionProblem(plasma, region, k), GRPF(; tol))
    # The genuine symmetric ±kvA root pair is recovered.
    @test count(ω -> abs(ω[]) > 0.02, sol.roots) ≥ 2
end
