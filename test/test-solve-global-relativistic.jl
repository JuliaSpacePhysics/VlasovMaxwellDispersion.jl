@testitem "GlobalDispersionProblem uses GRPF" begin
    plasma = Plasma(Species(0.0, 0.0, ColdVDF()))
    k = Wavenumber(0.0, 1.0)
    sol = solve(GlobalDispersionProblem(plasma, k, (0.4 - 0.2im, 1.4 + 0.2im)), GRPF(; tol=0.03))
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
    coldish = Species(-1.0, 0.2, Maxwellian(vth))
    rel = Species(-1.0, 0.2, MaxwellJuttner(mu))

    χcoldish = VM.contribution(coldish, omega, k)
    χrel = VM.contribution(rel, omega, k)

    @test maximum(abs.(χrel .- χcoldish)) / maximum(abs.(χcoldish)) < 1.0e-3
end

@testitem "Maxwell-Juttner trait and lower-half guard" begin
    using VlasovMaxwellDispersion
    using VlasovMaxwellDispersion: VlasovMaxwellDispersion as VM

    species = Species(-1.0, 0.2, MaxwellJuttner(1.0e4))

    @test VM.Regime(species) isa Relativistic
    @test all(isfinite, VM.contribution(Species(-1.0, 1.0, MaxwellJuttner(2.0)),
                                        1.0965 - 2.2732e-7im, Wavenumber(0.001, 0.1)))
    @test_throws ArgumentError VM.contribution(species, 0.8 - 0.01im, Wavenumber(0.1, 0.0))
end
