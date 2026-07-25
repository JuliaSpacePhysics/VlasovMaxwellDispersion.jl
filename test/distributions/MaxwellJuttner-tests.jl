@testitem "Maxwell-Juttner trait and lower-half guard" begin
    using VlasovMaxwellDispersion: VlasovMaxwellDispersion as VM

    species = NormalizedSpecies(-1.0, 0.2, MaxwellJuttner(1.0e4))

    @test VM.regime(species) isa Relativistic
    @test all(isfinite, contribution(species, 1.0965 - 2.2732e-7im, Wavenumber(0.001, 0.1)))
    @test_throws ArgumentError VM.contribution(species, 0.8 - 0.01im, Wavenumber(0.1, 0.0))
end

@testitem "Maxwell-Juttner tends to Maxwellian for large mu" begin
    let mu = 1.0e4, k = Wavenumber(0.1, 0.2), ω = 0.8 + 0.01im
        χrel = contribution(NormalizedSpecies(-1.0, 0.2, MaxwellJuttner(mu)), ω, k)
        χcoldish = contribution(NormalizedSpecies(-1.0, 0.2, Maxwellian(sqrt(2 / mu))), ω, k)
        @test χrel ≈ χcoldish rtol = 1.0e-3
    end
end
