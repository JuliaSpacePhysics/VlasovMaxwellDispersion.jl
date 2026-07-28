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

# Wrap each case-study page's @example execution in a hidden wall-clock timer.
const TIMINGS_FILE = joinpath(@__DIR__, "build-timings.tsv")
ENV["VMD_TIMINGS"] = TIMINGS_FILE
rm(TIMINGS_FILE; force = true)

timing_preprocess(name) = str ->
    "_t0 = time(); nothing #hide\n" * str *
    "\n#\nopen(ENV[\"VMD_TIMINGS\"], \"a\") do io; println(io, \"$name\\t\", round(time() - _t0; digits = 2)); end #hide\nnothing #hide\n"

foreach(LITERATE_SOURCES) do source
    name = first(splitext(basename(source)))
    Literate.markdown(source, CASE_STUDIES; documenter = true, preprocess = timing_preprocess(name))
end

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
    modules = [VlasovMaxwellDispersion, VlasovMaxwellDispersion.PlasmaBase, VlasovMaxwellDispersion.ReturnCode],
    authors = "Beforerr <zzj956959688@gmail.com> and contributors",
    sitename = "VlasovMaxwellDispersion.jl",
    format = Documenter.HTML(;
        canonical = "https://JuliaSpacePhysics.github.io/VlasovMaxwellDispersion.jl",
        size_threshold_ignore = HIDDEN_PAGES,          # embedded Typst report carries base64 figures
    ),
    checkdocs = :none,                                 # site is benchmark-focused, not full API ref
    pages = PAGES,
)

# Per-page build timing → CI step summary (or stdout locally), slowest first.
isfile(TIMINGS_FILE) && summary_timings(TIMINGS_FILE)

deploydocs(;
    repo = "github.com/JuliaSpacePhysics/VlasovMaxwellDispersion.jl",
    push_preview = true,
)
