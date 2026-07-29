# Report figures. Run after 01-04:
#   julia --project=docs experiments/eic-transverse-heating/05_make_figures.jl
using CairoMakie, DelimitedFiles, Printf

const DIR = @__DIR__
rd(n) = readdlm(joinpath(DIR, n), '\t'; skipstart = 1)
const TI = 0.3

# --- Fig 1: currentless growth, and the collapse of the current threshold ----
let ring = rd("out_ring.tsv"), joint = rd("out_joint.tsv")
    fig = Figure(size = (980, 400))

    ax1 = Axis(fig[1, 1]; xlabel = "O⁺ ring energy E⊥ [eV]", ylabel = "peak growth γ / Ω_O⁺",
        title = "no current anywhere: the ring alone (δ = 1)")
    E, g = ring[:, 2], ring[:, 3]
    lines!(ax1, E, max.(g, 0); color = :orangered, linewidth = 2.5)
    scatter!(ax1, E, max.(g, 0); color = :orangered, markersize = 7)
    vlines!(ax1, [0.63]; color = :gray, linestyle = :dash)
    text!(ax1, 0.70, 0.42; text = "threshold\nE⊥ = 0.63 eV", fontsize = 9, color = :gray)
    ylims!(ax1, -0.015, maximum(g) * 1.12)

    ax2 = Axis(fig[1, 2]; xlabel = "O⁺ ring energy E⊥ [eV]",
        ylabel = "threshold current j_c / j_c(unheated)",
        title = "what the heating does to the current requirement")
    j0 = 16.468
    hspan!(ax2, 1.0, 1.8; color = (:firebrick, 0.07))
    hspan!(ax2, -0.05, 1.0; color = (:seagreen, 0.07))
    hlines!(ax2, [1.0]; color = :black, linestyle = :dash, linewidth = 1.2)
    text!(ax2, 0.03, 1.68; text = "harder than with no heating", fontsize = 9, color = :firebrick)
    for (i, δ) in enumerate((0.1, 0.3, 1.0))
        m = joint[:, 1] .≈ δ
        Ed, r = joint[m, 3], joint[m, 6] ./ j0
        lines!(ax2, Ed, r; color = i, colormap = :viridis, colorrange = (1, 3.6),
            linewidth = 2.5, label = "ring fraction δ = $δ")
        scatter!(ax2, Ed, r; color = i, colormap = :viridis, colorrange = (1, 3.6), markersize = 8)
    end
    text!(ax2, 0.70, 0.07; text = "δ = 1: unstable at zero current", fontsize = 9, color = :black)
    axislegend(ax2; position = :rc, framevisible = false, labelsize = 9)
    ylims!(ax2, -0.06, 1.78)

    save(joinpath(DIR, "fig_collapse.png"), fig; px_per_unit = 2)
end

# --- Fig 2: the same moments, a different answer -----------------------------
let c = rd("out_control.tsv")
    fig = Figure(size = (980, 400))
    A, gr, gm, ur, um = c[:, 2], c[:, 3], c[:, 6], c[:, 5], c[:, 8]

    ax1 = Axis(fig[1, 1]; xlabel = "T⊥/T∥ (identical for both curves)",
        ylabel = "peak growth γ / Ω_O⁺ at zero current",
        title = "same n, T⊥, T∥ — different shape")
    lines!(ax1, A, max.(gr, 0); color = :orangered, linewidth = 2.5, label = "gyro-ring")
    scatter!(ax1, A, max.(gr, 0); color = :orangered, markersize = 8)
    lines!(ax1, A, max.(gm, 0); color = :steelblue, linewidth = 2.5, label = "bi-Maxwellian")
    scatter!(ax1, A, max.(gm, 0); color = :steelblue, markersize = 8)
    axislegend(ax1; position = :lt, framevisible = false, labelsize = 10)

    ax2 = Axis(fig[1, 2]; xlabel = "T⊥/T∥ (identical for both curves)",
        ylabel = "threshold drift u_c / v_e", title = "the moment-based threshold is not a threshold")
    lines!(ax2, A, ur; color = :orangered, linewidth = 2.5, label = "gyro-ring")
    scatter!(ax2, A, ur; color = :orangered, markersize = 8)
    lines!(ax2, A, um; color = :steelblue, linewidth = 2.5, label = "bi-Maxwellian")
    scatter!(ax2, A, um; color = :steelblue, markersize = 8)
    axislegend(ax2; position = :lb, framevisible = false, labelsize = 10)

    save(joinpath(DIR, "fig_moments.png"), fig; px_per_unit = 2)
end

println("wrote fig_collapse.png, fig_moments.png")
