using VlasovMaxwellDispersion
using Documenter
using Literate
using CairoMakie
using Typst_jll

DocMeta.setdocmeta!(VlasovMaxwellDispersion, :DocTestSetup, :(using VlasovMaxwellDispersion); recursive = true)

# Literate scripts → executed markdown pages (figures rendered at build time).
const SRC = joinpath(@__DIR__, "src")
const CASE_STUDIES = joinpath(SRC, "case-studies")
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

# Typst notes → Documenter pages
function typst_page(source, out_md; root = dirname(source))
    html = tempname() * ".html"
    run(`$(typst()) compile --features html --format html --root $root $source $html`)
    doc = read(html, String)
    title = match(r"<title>(.*?)</title>"s, doc).captures[1]
    style = match(r"<style>.*?</style>"s, doc).match
    body = match(r"<body>(.*)</body>"s, doc).captures[1]
    open(out_md, "w") do io
        println(io, "# ", title, "\n")
        println(io, "```@raw html\n", style, "\n", body, "\n```")
    end
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
const PAGES = vcat(ROOT_PAGES, ["Case studies" => CASE_STUDY_PAGES])

makedocs(;
    modules = [VlasovMaxwellDispersion, VlasovMaxwellDispersion.ReturnCode],
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
if isfile(TIMINGS_FILE)
    rows = sort!(split.(readlines(TIMINGS_FILE), '\t'); by = r -> -parse(Float64, r[2]))
    summary = get(ENV, "GITHUB_STEP_SUMMARY", nothing)
    io = isnothing(summary) ? stdout : open(summary, "a")
    println(io, "\n## Case-study build time\n\n| page | seconds |\n|---|---:|")
    foreach(r -> println(io, "| ", r[1], " | ", r[2], " |"), rows)
    isnothing(summary) || close(io)
end

deploydocs(;
    repo = "github.com/JuliaSpacePhysics/VlasovMaxwellDispersion.jl",
    push_preview = true,
)
