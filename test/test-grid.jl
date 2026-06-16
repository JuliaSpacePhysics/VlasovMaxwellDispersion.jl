@testitem "GridVDF(bi-Maxwellian) ≡ analytic Maxwellian (closed-form Z/Γ_n tensor)" begin
    using VlasovMaxwellDispersion: contribution
    vthp, vthq = 0.9, 1.2
    f0(u, v) = exp(-(u / vthp)^2) / (sqrt(pi) * vthp) * exp(-(v / vthq)^2) / (pi * vthq^2)
    vpar = collect(range(-6vthp, 6vthp, length=81))
    vperp = collect(range(0.0, 6vthq, length=61))
    F = [f0(u, v) for u in vpar, v in vperp]
    g = GridVDF(vpar, vperp, F; tol=1e-4)
    k = Wavenumber(0.1, 0.4)
    ω = 1.3 - 0.05im
    χg = contribution(Species(-1.0, 0.5, g), ω, k)
    χm = contribution(Species(-1.0, 0.5, Maxwellian(vth_par=vthp, vth_perp=vthq)), ω, k)
    acc3 = maximum(abs.(χg .- χm)) / maximum(abs.(χm))
    @test acc3 < 5e-3

    g = GridVDF(vpar, vperp, F; method=BicubicHermite(), tol=1e-4)
    χg = contribution(Species(-1.0, 0.5, g), ω, k)
    @test maximum(abs.(χg .- χm)) / maximum(abs.(χm)) < 6e-2

    # Non-cubic NonnegBSpline with order=4 ≥ cubic accuracy.
    g4 = GridVDF(vpar, vperp, F; method=NonnegBSpline{4}(tol=1e-4))
    χg = contribution(Species(-1.0, 0.5, g4), ω, k)
    acc4 = maximum(abs.(χg .- χm)) / maximum(abs.(χm))
    @test acc4 < 5e-3
    @test acc3 > acc4
end

# isotropic Maxwell–Jüttner sampled on a (p∥,p⊥) grid. (γmax clipped to the grid, support-zeroed bicubic, resonance-gated Plemelj)
@testitem "Relativistic GridVDF reproduces Maxwell–Jüttner" begin
    μ = 40.0
    γ(u, w) = sqrt(1 + u^2 + w^2)
    f0(u, w) = exp(-μ * γ(u, w))
    L = sqrt((1 + 25 / μ)^2 - 1)
    ppar = collect(range(-L, L, length=81))
    pperp = collect(range(0.0, L, length=61))
    F = [f0(u, w) for u in ppar, w in pperp]
    grel = Species(1.0, 1.0, GridVDF(ppar, pperp, F; tol=1e-4); regime=Relativistic())
    ref = Species(1.0, 1.0, MaxwellJuttner(mu=μ))
    k = Wavenumber(0.7, 0.4)
    for ω in (0.3 - 0.005im, 0.3 + 0.05im)
        χg = contribution(grel, ω, k)
        χr = contribution(ref, ω, k)
        @test maximum(abs.(χg .- χr)) / maximum(abs.(χr)) < 2e-3
    end
end

@testitem "GridVDF ≈ CoupledVDF on inseparable f₀" tags=[:slow] begin
    using VlasovMaxwellDispersion
    using VlasovMaxwellDispersion: contribution
    g0(u, v) = exp(-(u^2 + v^2 + 0.6u * v))

    L = 6.0
    vpar = collect(range(-L, L, length=81))
    vperp = collect(range(0.0, L, length=61))
    F = [g0(u, v) for u in vpar, v in vperp]
    g = GridVDF(vpar, vperp, F; tol=1e-4)
    cpl = CoupledVDF(g0; parlower=(-L), parupper=L, perpupper=L)
    k = Wavenumber(0.3, 0.4)
    ω = 1.2 - 0.05im
    χg = contribution(Species(-1.0, 1.0, g), ω, k)
    χc = contribution(Species(-1.0, 1.0, cpl), ω, k)
    @test maximum(abs.(χg .- χc)) / maximum(abs.(χc)) < 5e-3
end

# ref: test-coupled-external.jl
@testitem "GridVDF bi-kappa ≡ analytic bi-kappa" begin
    using VlasovMaxwellDispersion: contribution
    vA, κ = 1e-4, 6.0
    a2 = (2κ - 3) / (2κ)
    s = κ * a2 * vA^2
    f(u, v) = (1 + (u^2 + v^2) / s)^(-1 - κ)
    du(u, v) = (-1 - κ) * (1 + (u^2 + v^2) / s)^(-2 - κ) * (2u / s)
    dv(u, v) = (-1 - κ) * (1 + (u^2 + v^2) / s)^(-2 - κ) * (2v / s)
    L = 12vA
    vpar = collect(range(-L, L, length=81))
    vperp = collect(range(0.0, L, length=61))
    F = [f(u, v) for u in vpar, v in vperp]
    grid = Species(1.0, 1 / vA^2, GridVDF(vpar, vperp, F; tol=1e-4))
    exact = Species(1.0, 1 / vA^2, CoupledVDF(f; parlower=(-L), parupper=L, perpupper=L, dpar=du, dperp=dv))
    k = Wavenumber(1e-3 / vA, 0.03 / vA)         # k̃ = k_ALPS/vA = (300, 10)
    ω = 0.029311 - 9.9693e-6im                   # ALPS bi-kappa (κ=6) root_1
    χg = contribution(grid, ω, k)
    χe = contribution(exact, ω, k)
    @test maximum(abs.(χg .- χe)) / maximum(abs.(χe)) < 1e-2
end
