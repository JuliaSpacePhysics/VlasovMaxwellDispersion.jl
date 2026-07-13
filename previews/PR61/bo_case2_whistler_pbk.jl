# # Parallel whistler — electron product-bi-kappa (Guo et al. 2026, Case 2)
#
# Case 2 of the BO/ALPS solver comparison [arXiv:2606.14439](https://arxiv.org/abs/2606.14439):
# the parallel whistler instability driven by an anisotropic electron *product*
# bi-kappa `f_e ∝ [1 + v∥²/(κc∥²)]^{-(κ+1)}·[1 + v⊥²/(κc⊥²)]^{-(κ+1)}` with
# `T⟂e = 4 T∥e`, product-bi-kappa protons at `T_p = 50 eV`, for
# `κ ∈ {3, 7, ∞}` (paper Fig. 2). The `κ = 3` panel is where the paper reports
# BO's fitted growth rate deviating from the ALPS/benchmark curves.

using VlasovMaxwellDispersion
using DelimitedFiles, Printf
using CairoMakie

# ## Plasma setup
#
# SI parameters, normalized to the electron gyrofrequency `|ωce|` (electron-scale
# mode); velocities in units of `c`. The paper's temperature-preserving thermal
# speeds `c∥ = √(2qT∥/m·(1−1/2κ))`, `c⊥ = √(2qT⊥/m·(1−1/κ))` are exactly VMD's
# `ProductBiKappa` convention, so raw `vth = √(2qT/m)` maps both onto the same `f`.

const c0 = 2.99792458e8
qe = 1.602176634e-19; mp = 1.67262192369e-27; me = 9.1093837015e-31
eps0 = 8.8541878128e-12

B0 = 1.0e-8; n = 1.0e7
Tez = 102.0; Tep = 408.0; Tp = 50.0
κs = (3.0, 7.0, Inf)

wce = qe * B0 / me
Pi2e = n * qe^2 / (eps0 * me) / wce^2
Pi2p = n * qe^2 / (eps0 * mp) / wce^2
vthez = sqrt(2qe * Tez / me) / c0
vthep = sqrt(2qe * Tep / me) / c0
vthp = sqrt(2qe * Tp / mp) / c0

vdf(vth_para, vth_perp, κ) = isinf(κ) ?
    Maxwellian(; vth_para, vth_perp) :
    ProductBiKappa(; vth_para, vth_perp, kappa_para = κ, kappa_perp = κ)

# ## Seedless surveys, one per κ
#
# `k∥` is swept over `k·λₑ ∈ [0.01, 2.8]` (`λₑ = c/ωpe`, the paper's abscissa);
# the `ω` box spans the whistler band `0 < Re ω < |ωce|` up to the strongly
# growing peak `γ ≳ ωce` together with the damped branches around it.

kunit = sqrt(Pi2e)                     # k·λₑ → k c/|ωce|
region = (-0.2 - 0.2im, 1.0 + 1.5im)
geom = CartesianSweep(kz = vcat(0.01:0.005:0.3, 0.32:0.02:2.8) .* kunit)

sols = map(κs) do κ
    plasma = (
        NormalizedSpecies(me / mp, Pi2p, vdf(vthp, vthp, κ)),
        NormalizedSpecies(-1.0, Pi2e, vdf(vthez, vthep, κ)),
    )
    solve(GlobalDispersionProblem(plasma, region, geom))
end

# ## Verification against PlasmaBO
#
# The reference (`bo_case2_ref.tsv`) is the whistler branch tracked with
# PlasmaBO's *dedicated* product-bi-kappa solver `BOPBK` (`N = 2`) for
# `κ = 3, 7` — the low-κ-capable formulation of Bai et al. 2024, i.e. the
# paper's black benchmark curves, not the HH fit that fails at `κ = 3` — and
# with the Hermite–Hermite solver for the Maxwellian limit. At each reference
# `k` the nearest surveyed root is compared over the whole complex plane.

ref = readdlm(joinpath(@__DIR__, "bo_case2_ref.tsv"); comments = true)
kle(b) = [para(k) / kunit for k in b.k]
for (κ, sol) in zip(κs, sols)
    rows = ref[ref[:, 4] .== κ, :]
    Δmax = 0.0
    for r in eachrow(rows)
        r[3] > 0 || continue   # ref's damped tail exits the surveyed ω box
        ω_ref = complex(r[2], r[3])
        d = minimum(
            minimum((abs(ω - ω_ref) for (x, ω) in zip(kle(b), b.omega) if isfinite(ω) && abs(x - r[1]) < 0.005); init = Inf)
                for b in sol.roots
        )
        Δmax = max(Δmax, d)
    end
    @printf("κ=%-3s max |ω_vmd - ω_ref| = %.1e ωce\n", isinf(κ) ? "∞" : string(round(Int, κ)), Δmax)
end

# Both codes evaluate the analytic product-bi-kappa susceptibility for this
# parallel geometry, so the growing branch agrees to `~10⁻⁶ ωce` — far below
# the visual thickness of the paper's curves.
#
# !!! note "Convention discrepancy in the paper's κ = 3 panel"
#     With the paper's *written* thermal speeds `c∥ = √(2kT∥/m·(1−1/2κ))`,
#     `c⊥ = √(2kT⊥/m·(1−1/κ))` (used here and by PlasmaBO), the κ = 3 branch
#     crosses `γ = 0` at `k·λₑ ≈ 2.1`. The paper's Fig. 2(b) instead crosses at
#     `≈ 2.45`, which VMD reproduces only with *raw* `θ∥,⊥ = √(2kT/m)` (no κ
#     correction — evidently the convention of the plotted Bai 2024 benchmark).
#     At κ = 7 and ∞ the correction is within the curve width, so only the
#     κ = 3 panel differs visibly.

# ## Dispersion diagram — paper Fig. 2
#
# Left `Re ω`, right `γ`, rows `κ = 3, 7, ∞`; black dots: PlasmaBO track.
# The whistler branch grows for `k·λₑ ≈ 0.15–2` with peak `γ ≈ 1.08, 0.99,
# 0.92 |ωce|` at `k·λₑ ≈ 0.9` for `κ = 3, 7, ∞` — larger and at slightly larger
# `k` for smaller κ, matching the ALPS/benchmark curves (Fig. 2 b/d/f).

fig = Figure(size = (850, 780))
for (i, (κ, sol)) in enumerate(zip(κs, sols))
    lab = isinf(κ) ? "κ = ∞" : "κ = $(round(Int, κ))"
    axr = Axis(fig[i, 1]; ylabel = "Re ω / |ωce|", title = lab, xlabel = i == 3 ? "k λₑ" : "")
    axi = Axis(fig[i, 2]; ylabel = "γ / |ωce|", title = lab, xlabel = i == 3 ? "k λₑ" : "")
    for (j, b) in enumerate(sol.roots)
        x = kle(b)
        count(isfinite, b.omega) ≥ 3 || continue
        p = sortperm(x)
        lines!(axr, x[p], real.(b.omega)[p]; color = :royalblue, linewidth = 2)
        lines!(axi, x[p], imag.(b.omega)[p]; color = :orangered, linewidth = 2)
    end
    rows = ref[ref[:, 4] .== κ, :]
    scatter!(axr, rows[:, 1], rows[:, 2]; color = :black, markersize = 6)
    scatter!(axi, rows[:, 1], rows[:, 3]; color = :black, markersize = 6)
    hlines!(axi, [0.0]; color = (:black, 0.3), linestyle = :dash)
    xlims!(axr, 0, 2.6); xlims!(axi, 0, 2.6)
    ylims!(axr, 0, 1); ylims!(axi, -0.3, 1.5)
end
fig
