using VlasovMaxwellDispersion
using Documenter
using Literate
using CairoMakie

DocMeta.setdocmeta!(VlasovMaxwellDispersion, :DocTestSetup, :(using VlasovMaxwellDispersion); recursive = true)

# Literate scripts → executed markdown pages (figures rendered at build time).
const SRC = joinpath(@__DIR__, "src")
for f in ("cattaert.jl", "firehose_astfalk.jl", "ionbeam_gary84.jl")
    Literate.markdown(joinpath(SRC, f), SRC; documenter = true)
end

makedocs(;
    modules = [VlasovMaxwellDispersion],
    authors = "Beforerr <zzj956959688@gmail.com> and contributors",
    sitename = "VlasovMaxwellDispersion.jl",
    format = Documenter.HTML(;
        canonical = "https://JuliaSpacePhysics.github.io/VlasovMaxwellDispersion.jl",
    ),
    checkdocs = :none,                                 # site is benchmark-focused, not full API ref
    pages = [
        "Home" => "index.md",
        "Kappa-Maxwellian Plasma — Cattaert 2007" => "cattaert.md",
        "Firehose — bi-kappa (Astfalk 2017)" => "firehose_astfalk.md",
        "Ion beam (Gary 1984)" => "ionbeam_gary84.md",
    ],
)

deploydocs(;
    repo = "github.com/JuliaSpacePhysics/VlasovMaxwellDispersion.jl",
    push_preview = true,
)
