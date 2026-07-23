module VlasovMaxwellDispersionMakieExt

using VlasovMaxwellDispersion: dispersion_diagram!
import VlasovMaxwellDispersion
using Makie

# |k| for the x-axis — works for a k∥ sweep (m=1) or a (|k|,θ) sweep (m=2) alike.
_kmag(k) = sqrt(abs2(k))

struct FigureAxesPlots
    figure::Makie.Figure
    axes::Tuple{Makie.Axis,Makie.Axis}
    plots::Vector{Makie.AbstractPlot}
end
Base.show(io::IO, m::MIME, d::FigureAxesPlots) = show(io, m, d.figure)
Base.show(io::IO, m::MIME"text/plain", d::FigureAxesPlots) = show(io, m, d.figure)
Base.showable(m::MIME{M}, fg::FigureAxesPlots) where {M} = showable(m, fg.figure)
Base.display(d::FigureAxesPlots) = display(d.figure)
Base.iterate(d::FigureAxesPlots, s=1) = s > 3 ? nothing : (getfield(d, s), s + 1)

function VlasovMaxwellDispersion.dispersion_diagram!((axr, axi), sol; kw...)
    palette = Makie.wong_colors()
    plots = Makie.AbstractPlot[]
    for (i, root) in enumerate(sol.roots)
        col = palette[mod1(i, length(palette))]
        k = _kmag.(root.k)
        p = sortperm(k)
        k, ωr, ωi = k[p], real.(root.omega)[p], imag.(root.omega)[p]
        push!(plots, Makie.scatterlines!(axr, k, ωr; color=col, linewidth=2, markersize=4, kw...))
        push!(plots, Makie.scatterlines!(axi, k, ωi; color=col, linewidth=2, markersize=4, kw...))
    end
    return plots
end

"""
    dispersion_diagram(sol; figure=(;), title="")

Two axes side by side: `Re ω(k)` and `Im ω(k)`.

`fig, (axr, axi), plots = dispersion_diagram(sol)`.
"""
function VlasovMaxwellDispersion.dispersion_diagram(sol; figure=(;), xlabel="|k|", title="", kw...)
    fig = Makie.Figure(; figure...)
    axr = Makie.Axis(fig[1, 1]; xlabel, ylabel="Re ω")
    axi = Makie.Axis(fig[1, 2]; xlabel, ylabel="Im ω")
    isempty(title) || Makie.Label(fig[0, :], title; fontsize=16, font=:bold)
    Makie.linkxaxes!(axr, axi)
    plots = dispersion_diagram!((axr, axi), sol; kw...)
    Makie.colsize!(fig.layout, 1, Makie.Aspect(1, 1.0))
    Makie.colsize!(fig.layout, 2, Makie.Aspect(1, 1.0))
    Makie.resize_to_layout!(fig)
    return FigureAxesPlots(fig, (axr, axi), plots)
end

end
