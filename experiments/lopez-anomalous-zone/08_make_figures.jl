# Render report figures from figdata CSVs
using CairoMakie, DelimitedFiles

dir = @__DIR__
rd(n) = readdlm(joinpath(dir, n), ','; skipstart=1)

# --- Fig 1: holomorphy defect maps ---
hl, hc = rd("figdata_holo_lopez.csv"), rd("figdata_holo_corr.csv")
res, ims = sort(unique(hl[:, 1])), sort(unique(hl[:, 2]))
tomat(h) = [h[findfirst(r -> h[r, 1] == x && h[r, 2] == y, 1:size(h, 1)), 3] for x in res, y in ims]
fig = Figure(size=(880, 360))
for (i, (h, ttl)) in enumerate(((hl, "López Λ_L (with θ)"), (hc, "corrected continuation")))
    ax = Axis(fig[1, i]; xlabel="Re ω / Ωc", ylabel=i == 1 ? "Im ω / Ωc" : "", title=ttl)
    hm = heatmap!(ax, res, ims, tomat(h); colorrange=(-8, 0.5), colormap=:viridis)
    hlines!(ax, [0.0]; color=:white, linestyle=:dash, linewidth=1)
    scatter!(ax, [0.12415], [-0.13610]; marker=:xcross, color=:red, markersize=14)
    scatter!(ax, [0.16297], [-0.11048]; marker=:circle, color=:orangered, markersize=11)
    i == 2 && Colorbar(fig[1, 3], hm; label="log₁₀ |∂f/∂z̄| / |∂f/∂z|")
end
save(joinpath(dir, "fig_holomorphy.png"), fig)

# --- Fig 2: root trajectories, μ=2 ---
t = rd("figdata_traces.csv")
fig2 = Figure(size=(880, 360))
ax1 = Axis(fig2[1, 1]; xlabel="k∥ c / Ωc", ylabel="ωr / Ωc", title="root of each function vs k")
ax2 = Axis(fig2[1, 2]; xlabel="ωr / Ωc", ylabel="γ / Ωc", title="trajectory in the ω plane (k: 0.1→0.85)")
lines!(ax1, t[:, 1], t[:, 2]; color=:crimson, linewidth=2.5, label="López Λ_L (artifact descent)")
lines!(ax1, t[:, 1], t[:, 4]; color=:royalblue, linewidth=2.5, label="VMD det / corrected López")
axislegend(ax1; position=:lt, framevisible=false, labelsize=10)
lines!(ax2, t[:, 2], t[:, 3]; color=:crimson, linewidth=2.5)
lines!(ax2, t[:, 4], t[:, 5]; color=:royalblue, linewidth=2.5)
scatter!(ax2, [0.1630], [-0.1105]; marker=:diamond, color=:black, markersize=13, label="AAA pole (López UHP data)")
axislegend(ax2; position=:lb, framevisible=false, labelsize=10)
save(joinpath(dir, "fig_traces.png"), fig2)
println("figures written")

# Export figures from the public case study; run with --project=../../docs.
include(joinpath(@__DIR__, "../../docs/src/case-studies/relativistic_pair.jl"))
save(joinpath(@__DIR__, "fig_fig5_replica.png"), fig)