# In-band roots omitted from the public case study
using VlasovMaxwellDispersion
using CairoMakie

pair(vdf) = (NormalizedSpecies(1.0, 1.0, vdf), NormalizedSpecies(-1.0, 1.0, vdf))
p2 = pair(MaxwellJuttner(mu=2.0))
p10 = pair(MaxwellJuttner(mu=10.0))

k2 = collect(0.25:0.02:4.5)
ω2 = solve(DispersionProblem(p2, 0.13 - 0.227im, Wavenumber.(0.0, k2); mode=:L)).omega
k10a = collect(1.5:0.02:4.5)
ω10a = solve(DispersionProblem(p10, 1.482 - 0.2422im, Wavenumber.(0.0, k10a); mode=:L)).omega
k10b = collect(0.4:0.02:1.4)
ω10b = solve(DispersionProblem(p10, 0.1752 - 0.3852im, Wavenumber.(0.0, k10b); mode=:L)).omega

blu, red = Makie.wong_colors()[1], Makie.wong_colors()[6]
fig = Figure(size=(880, 380))
axr = Axis(fig[1, 1]; xlabel="k∥ c / |Ω|", ylabel="ωr / |Ω|", title="in-band roots")
axi = Axis(fig[1, 2]; xlabel="k∥ c / |Ω|", ylabel="γ / |Ω|", title="damping")
for (k, ω, color, label) in (
    (k2, ω2, blu, "μ = 2"),
    (k10a, ω10a, red, "μ = 10, member 1"),
    (k10b, ω10b, (red, 0.55), "μ = 10, member 2"),
)
    lines!(axr, k, real.(ω); color, linewidth=2, label)
    lines!(axi, k, imag.(ω); color, linewidth=2)
end
lines!(axr, 0:4, 0:4; color=(:black, 0.3), linestyle=:dash, label="ω = k∥")
axislegend(axr; position=:lt, framevisible=false, labelsize=9)
save(joinpath(@__DIR__, "fig_quasimodes.png"), fig)
