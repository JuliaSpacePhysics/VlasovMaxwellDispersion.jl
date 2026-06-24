@testitem "GlobalDispersionProblem uses GRPF" begin
    plasma = NormalizedSpecies(0.0, 0.0, ColdVDF())
    k = Wavenumber(0.0, 1.0)
    sol = solve(GlobalDispersionProblem(plasma, k, (0.4 - 0.2im, 1.4 + 0.2im)), GRPF(; tol = 0.03))
    roots, poles = sol.omega, sol.poles

    @test isempty(poles)
    @test any(abs.(roots .- (1 + 0im)) .< 0.05)
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
