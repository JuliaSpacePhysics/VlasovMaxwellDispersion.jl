@testitem "GlobalDispersionProblem uses GRPF" begin
    plasma = NormalizedSpecies(0.0, 0.0, ColdVDF())
    k = Wavenumber(0.0, 1.0)
    sol = solve(GlobalDispersionProblem(plasma, (0.4 - 0.2im, 1.4 + 0.2im), k), GRPF(; tol = 0.03))
    @test any(x -> abs(x[] - (1 + 0im)) < 0.05, sol.roots)
end
