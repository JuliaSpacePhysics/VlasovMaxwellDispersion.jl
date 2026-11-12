@testitem "Polishes |ω|≪1 roots (relative seed spread)" begin
    # Dense cold e-p plasma: EMIC/Alfvén root near ω≈8.6e-5 (≪ 1).
    mp_me = 1836.15
    plasma = (
        NormalizedPopulation(-1.0, 100.0, Cold()),
        NormalizedPopulation(1 / mp_me, 100.0 / mp_me, Cold()),
    )
    k = Wavevector(0.0, 0.04)
    sol = solve(DispersionProblem(plasma, 8.5e-5 + 0.0im, k), Muller())
    @test sol.retcode == ReturnCode.Success
    @test abs(sol.omega - 8.5e-5) < 1.0e-5   # stayed on low-ω branch
    @test sol.resid < 1.0e-8                 # polished to genuine root
end

@testitem "Polishes a whistler root next to the Ωp resonance, cold or warm" begin
    # just below Ωp, the L mode's cyclotron resonance is a pole of det 𝒟 that Muller's
    # quadratic model cannot cross
    mp_me = 1836.15
    ω = 0.00050945
    for vth in (0.0, 1.0e-3)
        vdf(v) = iszero(v) ? Cold() : BiMaxwellian(; vth = v)
        plasma = (
            NormalizedPopulation(-1.0, 100.0, vdf(vth)),
            NormalizedPopulation(1 / mp_me, 100.0 / mp_me, vdf(vth / sqrt(mp_me))),
        )
        sol = solve(DispersionProblem(plasma, 0.998ω + 0im, Wavevector(0.0, 0.157)), Muller())
        @test sol.retcode == ReturnCode.Success
        @test sol.omega ≈ ω rtol = 1.0e-4   # stayed on the whistler branch
        @test sol.resid < 1.0e-8
    end
end

@testitem "Seeded reduction: implicit det at one mode's seed, det between modes" begin
    using VlasovMaxwellDispersion: _seeded_mode, FullDet, ImplicitDet
    using VlasovMaxwellDispersion.StaticArrays: SMatrix
    using LinearAlgebra
    # Non-normal: least left/right singular vectors nearly orthogonal, as where modes coalesce.
    T = SMatrix{3, 3, ComplexF64}([1.0e-3 1 0; 0 2.0e-3 0; 0 0 1])
    f = _seeded_mode(FullDet(), T)
    @test f isa ImplicitDet
    @test abs(f(T)) ≈ last(svdvals(T)) rtol = 1.0e-6
    @test _seeded_mode(FullDet(), SMatrix{3, 3, ComplexF64}(Diagonal([1, 1.0e-3, 1.0e-3]))) === FullDet()
end
