@testitem "GridVDF(bi-Maxwellian) ≡ analytic Maxwellian (closed-form Z/Γ_n tensor)" begin
    vthp, vthq = 0.9, 1.2
    f0(u, v) = exp(-(u / vthp)^2) / (sqrt(pi) * vthp) * exp(-(v / vthq)^2) / (pi * vthq^2)
    vpar = range(-6vthp, 6vthp, length = 81)
    vperp = range(0.0, 6vthq, length = 61)
    F = [f0(u, v) for v in vperp, u in vpar]      # F[perp,para]
    g = GridVDF(vperp, vpar, F; rtol = 1.0e-4)
    k = Wavenumber(0.1, 0.4)
    ω = 1.3 - 0.05im
    χg = contribution(NormalizedSpecies(-1.0, 0.5, g), ω, k)
    χm = contribution(NormalizedSpecies(-1.0, 0.5, Maxwellian(vth_para = vthp, vth_perp = vthq)), ω, k)
    acc3 = maximum(abs.(χg .- χm)) / maximum(abs, χm)
    @test acc3 < 5.0e-3

    g = GridVDF(vperp, vpar, F; method = BicubicHermite())
    χg = contribution(NormalizedSpecies(-1.0, 0.5, g), ω, k)
    @test maximum(abs.(χg .- χm)) / maximum(abs, χm) < 6.0e-2

    # Non-cubic NonnegBSpline with order=4 ≥ cubic accuracy.
    g4 = GridVDF(vperp, vpar, F; method = NonnegBSpline{4}(rtol = 1.0e-4))
    χg = contribution(NormalizedSpecies(-1.0, 0.5, g4), ω, k)
    acc4 = maximum(abs.(χg .- χm)) / maximum(abs, χm)
    @test acc4 < 5.0e-3
    @test acc3 > acc4
end

# isotropic Maxwell–Jüttner sampled on a (p∥,p⊥) grid, routed through the coupled (p⊥,p∥) path
@testitem "Relativistic GridVDF reproduces Maxwell–Jüttner" begin
    μ = 40.0
    γ(u, w) = sqrt(1 + u^2 + w^2)
    f0(u, w) = exp(-μ * γ(u, w))
    L = sqrt((1 + 25 / μ)^2 - 1)
    ppar = range(-L, L, length = 81)
    pperp = range(0.0, L, length = 61)
    F = [f0(u, w) for w in pperp, u in ppar]      # F[perp,para]
    grel = GridVDF(pperp, ppar, F; rtol = 1.0e-4, regime = Relativistic())
    ref = MaxwellJuttner(mu = μ)
    k = Wavenumber(0.7, 0.4)
    for ω in (0.3 - 0.005im, 0.3 + 0.05im)
        χg = contribution(grel, ω, k)
        χr = contribution(ref, ω, k)
        @test χg ≈ χr rtol = 1.0e-3
    end
end

@testitem "GridVDF ≈ CoupledVDF on inseparable f₀" begin
    g0(u, v) = exp(-(u^2 + v^2 + 0.6u * v))

    L = 6.0
    vpar = range(-L, L, length = 81)
    vperp = range(0.0, L, length = 61)
    F = [g0(u, v) for v in vperp, u in vpar]      # F[perp,para]
    g = GridVDF(vperp, vpar, F; rtol = 1.0e-4)
    cpl = CoupledVDF(g0; para = (-L, L), perp = L)
    k = Wavenumber(0.3, 0.4)
    ω = 1.2 - 0.05im
    χc = contribution(NormalizedSpecies(-1.0, 1.0, cpl), ω, k)
    χg = contribution(NormalizedSpecies(-1.0, 1.0, g), ω, k)
    χg_cpl = contribution(NormalizedSpecies(-1.0, 1.0, g.coupled), ω, k)
    @test χg ≈ χc rtol = 5.0e-3
    @test χg_cpl ≈ χc rtol = 5.0e-3
end

# ref: test-coupled-external.jl
@testitem "GridVDF bi-kappa ≡ analytic bi-kappa" begin
    vA, κ = 1.0e-4, 6.0
    a2 = (2κ - 3) / (2κ)
    s = κ * a2 * vA^2
    f(w, u) = (1 + (u^2 + w^2) / s)^(-1 - κ)
    L = 12vA
    vpar = range(-L, L, length = 81)
    vperp = range(0.0, L, length = 61)
    F = [f(v, u) for v in vperp, u in vpar]       # F[perp,para]
    grid = NormalizedSpecies(1.0, 1 / vA^2, GridVDF(vperp, vpar, F; rtol = 1.0e-4))
    exact = NormalizedSpecies(1.0, 1 / vA^2, CoupledVDF(f; para = (-L, L), perp = L))
    k = Wavenumber(1.0e-3 / vA, 0.03 / vA)         # k̃ = k_ALPS/vA = (300, 10)
    ω = 0.029311 - 9.9693e-6im                   # ALPS bi-kappa (κ=6) root_1
    χg = contribution(grid, ω, k)
    χe = contribution(exact, ω, k)
    @test χg ≈ χe rtol = 1.0e-2
end
