using VlasovMaxwellDispersion
using TestItemRunner

@run_package_tests verbose=true
# @run_package_tests verbose=true filter=ti -> !(:slow in ti.tags)

@testitem "cold dispersion vs Stix" begin
    include("test-cold-stix.jl")
end
@testitem "SciML interop via residual seam" begin
    include("test-sciml-interop.jl")
end