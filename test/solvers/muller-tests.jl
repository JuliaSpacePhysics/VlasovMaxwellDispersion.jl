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
