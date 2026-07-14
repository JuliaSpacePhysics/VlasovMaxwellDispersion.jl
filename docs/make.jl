using VlasovMaxwellDispersion
using Documenter
using Literate
using CairoMakie

DocMeta.setdocmeta!(VlasovMaxwellDispersion, :DocTestSetup, :(using VlasovMaxwellDispersion); recursive = true)

# Literate scripts → executed markdown pages (figures rendered at build time).
const SRC = joinpath(@__DIR__, "src")
const CASE_STUDIES = joinpath(SRC, "case-studies")
const LITERATE_SOURCES = sort(filter(f -> endswith(f, ".jl"), readdir(CASE_STUDIES; join = true)))

foreach(LITERATE_SOURCES) do source
    Literate.markdown(source, CASE_STUDIES; documenter = true)
end

const ROOT_PAGES = sort(filter(f -> endswith(f, ".md"), readdir(SRC)); by = f -> (f != "index.md", f))
const CASE_STUDY_PAGES = map(
    f -> joinpath("case-studies", f),
    sort(filter(f -> endswith(f, ".md"), readdir(CASE_STUDIES))),
)
const PAGES = vcat(ROOT_PAGES, ["Case studies" => CASE_STUDY_PAGES])

makedocs(;
    modules = [VlasovMaxwellDispersion],
    authors = "Beforerr <zzj956959688@gmail.com> and contributors",
    sitename = "VlasovMaxwellDispersion.jl",
    format = Documenter.HTML(;
        canonical = "https://JuliaSpacePhysics.github.io/VlasovMaxwellDispersion.jl",
    ),
    checkdocs = :none,                                 # site is benchmark-focused, not full API ref
    pages = PAGES,
)

deploydocs(;
    repo = "github.com/JuliaSpacePhysics/VlasovMaxwellDispersion.jl",
    push_preview = true,
)
