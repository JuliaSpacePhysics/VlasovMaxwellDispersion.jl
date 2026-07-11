@testitem "GlobalDispersionProblem uses GRPF" begin
    plasma = NormalizedSpecies(0.0, 0.0, ColdVDF())
    k = Wavenumber(0.0, 1.0)
    sol = solve(GlobalDispersionProblem(plasma, (0.4 - 0.2im, 1.4 + 0.2im), k), GRPF(; tol = 0.03))
    @test any(x -> abs(x[] - (1 + 0im)) < 0.05, sol.roots)
end


@testitem "Maxwell-Juttner tends to Maxwellian for large mu" begin
    using VlasovMaxwellDispersion
    using VlasovMaxwellDispersion: VlasovMaxwellDispersion as VM

    mu = 1.0e4
    vth = sqrt(2 / mu)
    k = Wavenumber(0.1, 0.2)
    omega = 0.8 + 0.01im
    coldish = NormalizedSpecies(-1.0, 0.2, Maxwellian(vth))
    rel = NormalizedSpecies(-1.0, 0.2, MaxwellJuttner(mu))

    χcoldish = VM.contribution(coldish, omega, k)
    χrel = VM.contribution(rel, omega, k)

    @test maximum(abs.(χrel .- χcoldish)) / maximum(abs.(χcoldish)) < 1.0e-3
end
