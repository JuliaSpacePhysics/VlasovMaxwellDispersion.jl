# # Oblique proton firehose — coupled bi-kappa protons (Guo et al. 2026, Case 1)
#
# Case 1 of the BO/ALPS solver comparison [arXiv:2606.14439](https://arxiv.org/abs/2606.14439):
# the oblique proton firehose driven by an anisotropic *coupled* kappa distribution
# `f_p ∝ [1 + v∥²/(κc∥²) + v⊥²/(κc⊥²)]^{-(κ+1)}` with `T∥p = 2 T⟂p` at `θ = 45°`,
# Maxwellian electrons, for `κ ∈ {4, 6, 8}` (paper Fig. 1 b/d/f).

using VlasovMaxwellDispersion
using DelimitedFiles, Printf
using CairoMakie

# ## Plasma setup
#
# SI parameters, normalized to the proton gyrofrequency `ωcp`; velocities in
# units of `c`.

const c0 = 2.99792458e8
qe = 1.602176634e-19; mp = 1.67262192369e-27; me = 9.1093837015e-31
eps0 = 8.8541878128e-12; mu0 = 1.25663706212e-6

B0 = 0.1; n = 5.0e19
Te = 496.683; Tpz = 1986.734; Tpp = 993.367
κs = (4.0, 6.0, 8.0)

wcp = qe * B0 / mp
Pi2p = n * qe^2 / (eps0 * mp) / wcp^2
Pi2e = n * qe^2 / (eps0 * me) / wcp^2
vthpz = sqrt(2qe * Tpz / mp) / c0
vthpp = sqrt(2qe * Tpp / mp) / c0
vthe = sqrt(2qe * Te / me) / c0

# ## Seedless surveys, one per κ
#
# `k` is swept over `k·dᵢ ∈ [0.03, 0.45]` (`dᵢ = c/ωpp = vA/ωcp`, the paper's `λ_p`);
# the `ω` box straddles the real axis to capture the non-propagating firehose
# branch together with the damped modes around it.

vA = B0 / sqrt(mu0 * n * mp)
kunit = c0 / vA                        # k·dᵢ → k c/ωcp
region = (-0.1 - 0.25im, 0.5 + 0.12im)
geom = AngleSweep(k = collect(0.03:0.015:0.45) .* kunit, theta = deg2rad(45))

sols = map(κs) do κ
    vdf_p = LowRankVDF(
        BiKappa(vth_para = vthpz, vth_perp = vthpp, kappa = κ);
        rtol = 1.0e-10,
        para = (-10vthpz, 10vthpz),
        perp = 10vthpp,
    )
    plasma = (
        NormalizedSpecies(1.0, Pi2p, vdf_p),
        NormalizedSpecies(-mp / me, Pi2e, Maxwellian(vthe)),
    )
    solve(DispersionProblem(plasma, region, geom))
end

# ## Verification against PlasmaBO
#
# The reference (`bo_case1_ref.tsv`) is the unstable branch from PlasmaBO's
# Hermite–Hermite solver (`N = 2`, `J = 24`) for `κ = 6, 8`. At `κ = 4` the HH
# fit of the sampled kappa distribution is poor in the tails — the very effect
# the paper reports for BO (Fig. 1b) — so it is not a valid reference there and
# is omitted; VMD evaluates the analytic bi-kappa susceptibility and should
# instead track the paper's ALPS curve.

ref = readdlm(joinpath(@__DIR__, "bo_case1_ref.tsv"); comments = true)
kdi(b) = [sqrt(abs2(k)) / kunit for k in b.k]
## ref grid (Δk = 0.02) is offset from the survey grid (Δk = 0.015): nearest
## sample is ≤ 0.0075 away, so gate at half the survey spacing
γmax(sol, x0) = maximum(
    maximum((imag(ω) for (x, ω) in zip(kdi(b), b.omega) if isfinite(ω) && abs(x - x0) < 0.008); init = -Inf)
        for b in sol.roots
)
for (κ, sol) in zip(κs, sols)
    rows = ref[ref[:, 4] .== κ, :]
    isempty(rows) && continue
    Δmax = 0.0
    for r in eachrow(rows)
        r[3] > 0.005 || continue
        Δmax = max(Δmax, abs(γmax(sol, r[1]) - r[3]))
    end
    @printf("κ=%.0f  max |γ_vmd - γ_ref| = %.1e ωcp\n", κ, Δmax)
end

# Agreement at the few-`10⁻³ ωcp` truncation level of the reference's
# Hermite–Hermite expansion, as for the Astfalk case.

# ## Growth rates — paper Fig. 1 (b), (d), (f)
#
# Colored: all surveyed branches; black dots: PlasmaBO track (`κ = 6, 8`).
# Peak growth `γ ≈ 0.05–0.066 ωcp` near `k·dᵢ ≈ 0.25` grows with `κ`,
# matching the ALPS (blue dotted) curves of the paper.

fig = Figure(size = (700, 780))
palette = Makie.wong_colors()
for (i, (κ, sol)) in enumerate(zip(κs, sols))
    ax = Axis(
        fig[i, 1]; ylabel = "γ / ωcp", title = "κ = $(round(Int, κ))",
        xlabel = i == 3 ? "k dᵢ" : ""
    )
    for (j, b) in enumerate(sol.roots)
        x = kdi(b)
        p = sortperm(x)
        lines!(ax, x[p], imag.(b.omega)[p]; color = palette[mod1(j, length(palette))], linewidth = 2)
    end
    rows = ref[ref[:, 4] .== κ, :]
    isempty(rows) || scatter!(ax, rows[:, 1], rows[:, 3]; color = :black, markersize = 6)
    hlines!(ax, [0.0]; color = (:black, 0.3), linestyle = :dash)
    xlims!(ax, 0, 0.6)
    ylims!(ax, -0.02, 0.075)
end
fig
