# # Kappa-Maxwellian Plasma — Cattaert 2007
#
# A single electron population with strong anisotropy along the field
# (`κ∥ = 1`) and nearly Maxwellian across it (`κ⟂ = 200`).
# We first track four known electromagnetic branches from a single seed,
# then the seedless global solver rediscovers them all at once.
# Verified against [PlasmaBO.jl `rlp_Cattaert07` case](https://juliaspacephysics.github.io/PlasmaBO.jl/dev/rlp_Cattaert07/).

using VlasovMaxwellDispersion
using DelimitedFiles, Printf, Unitful

# ## Plasma setup
# Here we specify the temperature `T` and resolve `vth` from particle identity `Electron()` and `T`.
T = 2555 * u"eV" # vth = √(2qT/mₑ) = 0.1 c
plasma = Plasma(
    Species(Electron(), sc -> ProductBiKappa(sc; kappa_para=1, kappa_perp=200.0); n=2.43e6u"m^-3", T),
    B0=1.0e-6u"T",
)

# ## Reference roots and wavevectors
#
# The tabulated BOPBK roots (`cattaert07_ref.tsv`) use reference `kρ = k·vtp/ωce`;
# we invert this to VMD's `k` (normalized to `ωce/c`) at the fixed propagation angle.
# `vtp` is the perpendicular kappa `θ⊥ = √(1−1/κ⊥)·vth`.

vtp = sqrt(1 - 1 / 200) * scales(plasma).species[1].vth
ref = readdlm(joinpath(@__DIR__, "cattaert07_ref.tsv"); comments=true)
kρs = unique(ref[:, 1])
θ = deg2rad(30)
ks = [Wavenumber(kρ / vtp * sin(θ), kρ / vtp * cos(θ)) for kρ in kρs]
i0 = argmin(abs.(kρs .- 0.1))          # seed from kρ ≈ 0.1

# ## Seeded tracking of the four branches
#
# Each branch is followed bidirectionally from its `kρ ≈ 0.1` seed, then compared
# point-by-point against the BOPBK reference.

ttrack = @elapsed ωs = map(1:4) do ib
    rows = ref[ref[:, 2] .== ib, :]
    seed = complex(rows[i0, 3], rows[i0, 4])
    ω = solve(DispersionProblem(plasma, Seed(seed, ks[i0]), ks)).omega  # fans out both ways from kρ≈0.1
    ωref = complex.(rows[:, 3], rows[:, 4])
    dre = abs.(abs.(real.(ω)) .- abs.(real.(ωref)))
    dim = abs.(imag.(ω) .- imag.(ωref))
    @printf(
        "branch %d: seed=%.5f%+.2eim  maxΔ|Re|=%.2e  maxΔIm=%.2e  nfinite=%d/%d\n",
        ib, real(seed), imag(seed), maximum(dre), maximum(dim), count(isfinite, ω), length(ω)
    )
    ω
end
@printf("seeded tracking: %.1f s for 4 branches × %d k-points\n", ttrack, length(kρs))

# VMD reproduces the reference roots, including the weak damping.
#
# Solid lines: VMD tracks. Open markers: BOPBK reference.

using CairoMakie
fig = Figure(size=(700, 620))
axr = Axis(fig[1, 1]; ylabel="Re ω / ωce", title="Cattaert 2007 — VMD vs BOPBK")
axi = Axis(fig[2, 1]; xlabel="k ρce", ylabel="Im ω / ωce")
palette = Makie.wong_colors()
for ib in 1:4
    rows = ref[ref[:, 2] .== ib, :]
    col = palette[ib]
    lines!(axr, kρs, real.(ωs[ib]); color=col, linewidth=2, label="branch $ib")
    lines!(axi, kρs, imag.(ωs[ib]); color=col, linewidth=2)
    scatter!(axr, rows[:, 1], rows[:, 3]; color=col, markersize=5, marker=:circle)
    scatter!(axi, rows[:, 1], rows[:, 4]; color=col, markersize=5, marker=:circle)
end
ylims!(axr, 0, 3)
axislegend(axr; position=:lt, framevisible=false)
fig

# ## Seedless survey
#
# The same branches, discovered *without* initial points: a
# `DispersionProblem` over an `ω` box and a `kρ` scan finds all roots of
# `det 𝒟 = 0` at once. `region` is the `ω` search box (in units of `ωce`);
# `geom` sweeps `k` at fixed `θ`, spanning `kρ ∈ [0.005, 0.3]`.

region = (0.005 - 0.16im, 3.05 + 0.02im)
geom = AngleSweep(k=(0.05 / vtp * 0.1, 0.3 / vtp), theta=θ)
prob = DispersionProblem(plasma, region, geom)
sol = solve(prob)

# `dispersion_diagram` plots the surveyed branches: `Re ω(kρ)` and `Im ω(kρ)`,
# one colour per discovered branch. The four tabulated branches appear as
# continuous curves spanning the sweep; the extra flat branches near `ω ≈ ωce,
# 2ωce` are genuine but heavily damped kinetic roots (`Im ω ≲ −0.03`)
# Black dots overlay the BOPBK reference to confirm the discovered branches sit
# on the tabulated roots; the diagram's x-axis is `|k|` in units of `ωce/c`, so
# the reference `kρ` maps to `kρ/vtp`.

figs, (axr2, axi2), = dispersion_diagram(sol; title="Cattaert 2007 — seedless m=1 survey")
scatter!(axr2, ref[:, 1] ./ vtp, ref[:, 3]; color=:black, markersize=4)
scatter!(axi2, ref[:, 1] ./ vtp, ref[:, 4]; color=:black, markersize=4)
ylims!(axi2, -0.02, 0.005)
figs
