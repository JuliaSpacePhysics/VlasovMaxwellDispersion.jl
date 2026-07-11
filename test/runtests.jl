using VlasovMaxwellDispersion
using TestItemRunner

@run_package_tests verbose = true

@testitem "besselj_ladder vs SpecialFunctions" begin
    include("test-bessel-ladder.jl")
end
@testitem "cold dispersion vs Stix" begin
    include("test-cold-stix.jl")
end
@testitem "SciML interop" begin
    include("test-sciml-interop.jl")
end

@testitem "Aqua" begin
    using Aqua
    using VlasovMaxwellDispersion

    Aqua.test_all(VlasovMaxwellDispersion)
end
