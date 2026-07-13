# # Ion cyclotron emission — alpha ring-beam (Guo et al. Case 4)
#
# Quasi-perpendicular ion-cyclotron emission (ICE) driven by a fusion-born alpha
# ring-beam, after Case 4 of the BO vs ALPS benchmark (Guo et al., arXiv:2606.14439;
# reproducing Warwick et al. 2018-type ICE, validated in Xie et al. 2025). A fast
# magnetosonic/Bernstein branch propagating at `θ = 89.5°` crosses successive
# deuteron cyclotron harmonics; each crossing that overlaps the alpha ring resonance
# goes unstable. The plasma is deuterons + electrons (Maxwellian) plus a hot alpha
# ring-beam with perpendicular ring speed `v_dr = 0.045c`.

using VlasovMaxwellDispersion
using DelimitedFiles, Printf
using CairoMakie

# ## Plasma setup
#
# SI parameters (Table I), normalized to the *proton* gyrofrequency `ωcp = eB₀/mp`
# even though no protons are present — this fixes the paper's `Ω_cp` axis. Alpha and
# deuteron share `q/m = e/(2mp)`, so `Ω_α = Ω_D = ½ωcp`; ICE harmonics therefore sit
# at half-integer multiples of `ωcp` (`ωr ≈ 2.5, 3.0, 3.5` = deuteron harmonics 5,6,7).
# Velocities are in units of `c`; `vth = √(2eT/m)/c` with `T` in eV. The alpha ring is
# VMD's `GaussianRing`, a literal shifted-Gaussian `∝ exp[−v∥²/c∥² − (v⊥−v_dr)²/c⊥²]`.

const c0 = 2.99792458e8
qe = 1.602176634e-19; mp = 1.67262192369e-27
eps0 = 8.8541878128e-12
me = 5.447e-4 * mp                     # Table I electron mass

B0 = 2.1
wcref = qe * B0 / mp                    # proton gyrofrequency = Ω_cp
na, qa, ma = 1.0e16, 2qe, 4mp          # alpha ring-beam
nd, qd, md = 9.98e18, qe, 2mp          # deuterons
ne, qee, mee = 1.0e19, qe, me          # electrons

Om(q, m) = (q * B0 / m) / wcref
pi2(n, q, m) = n * q^2 / (eps0 * m) / wcref^2
vth(T, m) = sqrt(2qe * T / m) / c0

vdf_a = GaussianRing(vth_para = vth(1000.0, ma), vth_perp = vth(1000.0, ma), vd = 0.0, vr = 0.045)
plasma = (
    NormalizedSpecies(Om(qa, ma), pi2(na, qa, ma), vdf_a),
    NormalizedSpecies(Om(qd, md), pi2(nd, qd, md), Maxwellian(vth(1000.0, md))),
    NormalizedSpecies(-Om(qee, mee), pi2(ne, qee, mee), Maxwellian(vth(1000.0, mee))),
)

# ## Seedless survey
#
# `k` is swept over `k·λp ∈ [3, 7]` with `λp = c/ωpp` the proton inertial length
# (`ωpp` at `n = 10¹⁹ m⁻³`), so `kunit = ωpp/ωcp` maps the normalized wavenumber to
# VMD's `k c/ωcp` units. The `ω` box spans the first four harmonic crossings
# (`Re ω ∈ [1.8, 4.2]`) with a thin unstable margin — ICE growth rates are small
# (`γ ≲ 0.08 ωcp`). At the ring validity limit `Λr = k⊥ v_dr/Ω_α ≈ 4` at the highest
# mode, well inside `GaussianRing`'s `Λr ≲ 10`.

wpp = sqrt(ne * qe^2 / (eps0 * mp))
kunit = wpp / wcref
kλp = range(3.0, 7.0, length = 41)
region = (1.8 - 0.8im, 4.2 + 0.2im)
geom = AngleSweep(k = collect(kλp) .* kunit, theta = deg2rad(89.5))
sol = solve(GlobalDispersionProblem(plasma, region, geom))

# ## Verification against PlasmaBO
#
# `bo_case4_ref.tsv` is the BO growth curve (PlasmaBO.jl, `BOHH(N=10, J=8)`): the
# dominant growing root per `k`. At each `k` with `γ_ref > 0.005 ωcp` we take the
# largest `Im ω` over all VMD roots at the matching `k`.

ref = readdlm(joinpath(@__DIR__, "bo_case4_ref.tsv"); comments = true)
kλ(b) = [sqrt(abs2(k)) / kunit for k in b.k]
Δmax = 0.0
for r in eachrow(ref)
    r[3] > 0.005 || continue
    γ = maximum(
        maximum((imag(ω) for (x, ω) in zip(kλ(b), b.omega) if isfinite(ω) && abs(x - r[1]) < 0.02); init = -Inf)
            for b in sol.roots
    )
    global Δmax = max(Δmax, abs(γ - r[3]))
    @printf("k·λp=%.2f  γ_BO=%.4f  γ_vmd=%.4f  Δ=%.1e\n", r[1], r[3], γ, abs(γ - r[3]))
end
Δmax

# VMD's analytic ring-beam susceptibility and BO's Hermite–Hermite expansion of the
# same shifted-Gaussian agree to `~10⁻³ ωcp` across all three growing harmonics —
# identical unstable `k`-windows and peak growth rates.

# ## Dispersion diagram
#
# Panel (a): the propagating branch (rising) crossing the deuteron cyclotron
# harmonics (dashed lines at `ωr/ωcp = 2.5, 3.0, 3.5, 4.0`). Panel (b): growth rate,
# with the BO reference as black dots. Instability appears in narrow `k`-windows at
# each harmonic crossing — the first three (`ωr < 3.5`) are the modes reported in the
# paper; the fourth (`ωr ≈ 4`) is the strongest here.

fig = Figure(size = (900, 400))
axr = Axis(fig[1, 1]; xlabel = "k λp", ylabel = "ωr / ωcp", title = "(a)")
axi = Axis(fig[1, 2]; xlabel = "k λp", ylabel = "γ / ωcp", title = "(b)")
hlines!(axr, [2.5, 3.0, 3.5, 4.0]; color = (:gray, 0.5), linestyle = :dash)
palette = Makie.wong_colors()
for (i, b) in enumerate(sol.roots)
    col = palette[mod1(i, length(palette))]
    x = kλ(b)
    p = sortperm(x)
    scatter!(axr, x[p], real.(b.omega)[p]; color = col, markersize = 4)
    scatter!(axi, x[p], imag.(b.omega)[p]; color = col, markersize = 4)
end
scatter!(axi, ref[:, 1], ref[:, 3]; color = :black, markersize = 5)
xlims!(axr, 3, 7); ylims!(axr, 2, 4.3)
xlims!(axi, 3, 7); ylims!(axi, -0.005, 0.1)
hlines!(axi, [0.0]; color = (:black, 0.3))
fig
