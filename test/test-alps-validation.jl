@testitem "ALPS test_kpar_fast scan" begin
    using VlasovMaxwellDispersion
    using DelimitedFiles

    vA = 1.0e-4
    me = 5.44662e-4
    ion = NormalizedSpecies(1.0, 1 / vA^2, Maxwellian(vA))
    electron = NormalizedSpecies(-1 / me, 1 / (me * vA^2), Maxwellian(vA / sqrt(me)))
    plasma = (ion, electron)

    data = readdlm(joinpath(@__DIR__, "fixtures/alps/test_kpar_fast.scan_kpara_1.root_1"))
    ks = [Wavenumber(row[1] / vA, row[2] / vA) for row in eachrow(data)]
    reference = complex.(data[:, 3], data[:, 4])

    roots = solve(BranchProblem(plasma, ks, reference[1]), ArcLength(; atol = 1.0e-9, maxiter = 100)).omega

    @test all(isfinite, roots)
    @test maximum(abs.(real.(roots) .- real.(reference)) ./ abs.(real.(reference))) < 1.0e-2
    @test maximum(abs.(imag.(roots) .- imag.(reference)) ./ max.(abs.(imag.(reference)), 1.0e-12)) < 5.0e-2
end

@testitem "ALPS test_relativistic generated root 2" begin
    using VlasovMaxwellDispersion
    using DelimitedFiles
    using LinearAlgebra

    ion = NormalizedSpecies(1.0, 1.0, MaxwellJuttner(2.0))
    electron = NormalizedSpecies(-1.0, 1.0, MaxwellJuttner(2.0))
    plasma = (ion, electron)

    data = readdlm(joinpath(@__DIR__, "fixtures/alps/test_relativistic.scan_kpara_1.root_2"))
    residuals = map(eachrow(data)) do row
        k = Wavenumber(row[1], row[2])
        omega = complex(row[3], row[4])
        abs(det(𝒟(plasma, omega, k)))
    end

    @test maximum(residuals) < 1.0e-5
end

@testitem "ALPS test_relativistic root 1 (low branch)" begin
    using VlasovMaxwellDispersion
    using DelimitedFiles

    ion = NormalizedSpecies(1.0, 1.0, MaxwellJuttner(2.0))
    electron = NormalizedSpecies(-1.0, 1.0, MaxwellJuttner(2.0))
    plasma = (ion, electron)

    data = readdlm(joinpath(@__DIR__, "fixtures/alps/test_relativistic.scan_kpara_1.root_1"))

    # Raw |det D| is inflated ~1/ω⁴ by the curl-curl term at low ω̃ branch, so
    # confirm a genuine root via the scale-invariant sol.resid.
    for row in eachrow(data)
        k = Wavenumber(row[1], row[2])
        ωref = complex(row[3], row[4])
        sol = solve(LocalDispersionProblem(plasma, k, ωref))
        ω = sol.omega
        @test sol.resid < 1.0e-6
        # Exact analytic Maxwell-Jüttner vs ALPS's coarse 60×30 gridded+order-30
        # Damping (Im ω̃) of this near-marginal mode is far more approximation-sensitive
        # so only Re is asserted.
        @test real(ω) ≈ real(ωref) rtol = 2.0e-2
    end
end
