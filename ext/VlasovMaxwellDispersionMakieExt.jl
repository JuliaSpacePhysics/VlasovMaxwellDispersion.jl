module VlasovMaxwellDispersionMakieExt

using VlasovMaxwellDispersion: SurveySolution
import VlasovMaxwellDispersion
using Makie

# |k| for the x-axis — works for a k∥ sweep (m=1) or a (|k|,θ) sweep (m=2) alike.
_kmag(k) = sqrt(abs2(k))

"""
    dispersion_diagram(sol; figure=(;), title="dispersion survey") -> Figure

Two stacked axes: `Re ω(|k|)` and `Im ω(|k|)`. Root branches are drawn as lines
(sorted by `|k|`), pole branches as scatter. Each branch gets its own colour.
"""
function VlasovMaxwellDispersion.dispersion_diagram(
        sol::SurveySolution;
        figure = (;), title = "dispersion survey", kw...
    )
    fig = Makie.Figure(; figure...)
    axr = Makie.Axis(fig[1, 1]; ylabel = "Re ω", title)
    axi = Makie.Axis(fig[2, 1]; xlabel = "|k|", ylabel = "Im ω")
    Makie.linkxaxes!(axr, axi)
    palette = Makie.wong_colors()
    for (i, root) in enumerate(sol.roots)
        col = palette[mod1(i, length(palette))]
        k = _kmag.(root.k)
        p = sortperm(k)
        k, ωr, ωi = k[p], real.(root.omega)[p], imag.(root.omega)[p]
        Makie.scatterlines!(axr, k, ωr; color = col, linewidth = 2, markersize = 4, kw...)
        Makie.scatterlines!(axi, k, ωi; color = col, linewidth = 2, markersize = 4, kw...)
    end
    return fig
end

end
