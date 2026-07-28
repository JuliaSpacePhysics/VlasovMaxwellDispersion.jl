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

function summary_timings(file)
    rows = sort!(split.(readlines(file), '\t'); by = r -> -parse(Float64, r[2]))
    summary = get(ENV, "GITHUB_STEP_SUMMARY", nothing)
    io = isnothing(summary) ? stdout : open(summary, "a")
    println(io, "\n## Build time\n\n| page | seconds |\n|---|---:|")
    foreach(r -> println(io, "| ", r[1], " | ", r[2], " |"), rows)
    isnothing(summary) || close(io)
end