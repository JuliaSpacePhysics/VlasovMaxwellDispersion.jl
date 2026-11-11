using VlasovMaxwellDispersion
using Documenter
using Literate
using CairoMakie
using Typst_jll
include("utils.jl")

DocMeta.setdocmeta!(VlasovMaxwellDispersion, :DocTestSetup, :(using VlasovMaxwellDispersion); recursive = true)

# Literate scripts → executed markdown pages (figures rendered at build time).
const SRC = joinpath(@__DIR__, "src")
const CASE_STUDIES = joinpath(SRC, "case-studies")
const REPRESENTATIONS = joinpath(SRC, "representations")
const LITERATE_SOURCES = sort(filter(f -> endswith(f, ".jl"), readdir(CASE_STUDIES; join = true)))

foreach(source -> Literate.markdown(source, CASE_STUDIES; documenter = true), LITERATE_SOURCES)

const TYP_PAGES = [
    joinpath(SRC, "relativistic.typ") => joinpath(SRC, "relativistic.md"),
    joinpath(@__DIR__, "..", "experiments", "lopez-anomalous-zone", "report.typ") =>
        joinpath(SRC, "lopez-anomalous-zone.md"),
]
foreach(((src, md),) -> typst_page(src, md), TYP_PAGES)

const HIDDEN_PAGES = ["lopez-anomalous-zone.md"]
const ROOT_PAGES = sort(
    filter(f -> endswith(f, ".md") && f ∉ HIDDEN_PAGES, readdir(SRC));
    by = f -> (f != "index.md", f),
)
const CASE_STUDY_PAGES = map(
    f -> joinpath("case-studies", f),
    sort(filter(f -> endswith(f, ".md"), readdir(CASE_STUDIES))),
)
const REPRESENTATION_PAGES = map(
    f -> joinpath("representations", f),
    sort(filter(f -> endswith(f, ".md"), readdir(REPRESENTATIONS))),
)
const PAGES = vcat(
    ROOT_PAGES,
    ["VDF representations" => REPRESENTATION_PAGES],
    ["Case studies" => CASE_STUDY_PAGES],
)

makedocs(;
    modules = [VlasovMaxwellDispersion],
    authors = "Beforerr <zzj956959688@gmail.com> and contributors",
    sitename = "VlasovMaxwellDispersion.jl",
    format = Documenter.HTML(;
        canonical = "https://JuliaPlasma.github.io/VlasovMaxwellDispersion.jl",
        size_threshold_ignore = HIDDEN_PAGES,          # embedded Typst report carries base64 figures
    ),
    checkdocs = :none,                                 # site is benchmark-focused, not full API ref
    pages = PAGES,
)

deploydocs(;
    repo = "github.com/JuliaPlasma/VlasovMaxwellDispersion.jl",
    push_preview = true,
)
