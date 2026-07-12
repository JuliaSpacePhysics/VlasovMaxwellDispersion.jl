# Relativistic path: an isotropic Maxwell–Jüttner f₀ fed through the
# general CoupledVDF must reproduce the closed Maxwell–Jüttner (Swanson) tensor
@testitem "Relativistic CoupledVDF reproduces Maxwell–Jüttner" begin
    μ = 40.0
    L = sqrt((1 + 25 / μ)^2 - 1)
    ref = MaxwellJuttner(mu = μ)
    rel = CoupledVDF(ref; para = (-L, L), perp = L, regime = Relativistic())

    for ω in (0.3 - 0.005im, 0.3 + 0.05im), kperp in (0.0, 0.3, 0.6, 0.7)
        k = Wavenumber(kperp, 0.4)
        χA = contribution(rel, ω, k; closure = Newberger())
        χB = contribution(rel, ω, k)
        χref = contribution(ref, ω, k)
        @test χA ≈ χref rtol = 1.0e-5
        @test χB ≈ χref rtol = 1.0e-5
    end

    for ω in (0.3 - 0.05im), kperp in (0.0, 0.3, 0.6, 0.7)
        k = Wavenumber(kperp, 0.4)
        χA = contribution(rel, ω, k; closure = Newberger())
        χB = contribution(rel, ω, k)
        χref = contribution(ref, ω, k)
        @test χA ≈ χref rtol = 1.0e-5
        @test χB ≈ χref
    end
end

@testitem "Relativistic CoupledVDF finite at large k⊥ (off-support poles)" begin
    μ = 40.0
    γ(w, u) = sqrt(1 + u^2 + w^2)
    f0(w, u) = exp(-μ * γ(w, u))
    L = sqrt((1 + 25 / μ)^2 - 1)
    kw = (para = (-L, L), perp = L)
    s = CoupledVDF(f0; kw..., regime = Relativistic())
    ref = MaxwellJuttner(mu = μ)
    ω = 0.3 - 0.05im
    for kperp in (1.2, 2.0, 3.5)
        k = Wavenumber(kperp, 0.4)
        χB = contribution(s, ω, k)
        χA = contribution(s, ω, k; closure = Newberger())
        χref = contribution(ref, ω, k)
        @test χA ≈ χref rtol = 1.0e-4
        @test χB ≈ χref
    end
end

# Strongly-relativistic (μ=2) regime with resonances INSIDE the support
@testitem "Relativistic CoupledVDF matches Swanson at μ=2 (subluminal, any damping)" begin
    μ = 2.0
    ref = MaxwellJuttner(mu = μ)
    P = 8.0
    rel = CoupledVDF(ref; para = (-P, P), perp = P, regime = Relativistic())

    pperp = range(0.0, 5.0, length = 61)
    ppar = range(-5.0, 5.0, length = 121)
    F = [ref(w, u) for w in pperp, u in ppar]
    grel = GridVDF(pperp, ppar, F; rtol = 1.0e-4, regime = Relativistic())

    k = Wavenumber(0.7, 0.4)
    for ω in (0.3 + 0.0im, 0.3 + 0.05im, 0.3 - 0.001im, 0.3 - 0.005im, 0.3 - 0.02im, 0.3 - 0.1im, 0.3 - 0.15im)
        χB = contribution(rel, ω, k)
        χA = contribution(rel, ω, k; closure = Newberger())
        χg = contribution(grel, ω, k)
        χref = contribution(ref, ω, k)
        @test χB ≈ χref rtol = 2.0e-4
        @test χA ≈ χref rtol = 2.0e-4
        @test_broken χg ≈ χref rtol = 1.0e-3
    end
    # parallel propagation and oblique, damped
    for (kk, ω) in ((Wavenumber(0.0, 0.4), 0.3 - 0.01im), (Wavenumber(0.2, 0.6), 0.5 - 0.02im))
        χB = contribution(rel, ω, kk)
        χref = contribution(ref, ω, kk)
        @test χB ≈ χref rtol = 1.0e-4
    end
end

# Superluminal (|Re ω|>|k∥|): for Im ω ≥ 0 the straight (p⊥,p∥) integral is exact —
# including exactly-real ω, where the resonance-ellipse apex is an integrable kink.
# For Im ω < 0 the apex branch point crosses the p⊥ path: unsupported (warns);
# the way to reach damped superluminal ω is external continuation from Im ω ≥ 0
# samples — the recipe below (two sample grids must agree on the extrapolation).
@testitem "Relativistic superluminal: exact for Im ω ≥ 0, warns below" begin
    μ = 6.0
    P = 5.0
    ref = MaxwellJuttner(mu = μ)
    rel = NormalizedSpecies(1.0, 1.0, CoupledVDF(ref; para = (-P, P), perp = P, regime = Relativistic()))
    refs = NormalizedSpecies(1.0, 1.0, ref)
    k = Wavenumber(1.5, 0.8)
    for ω in (0.9 + 0.05im, 0.9 + 0.0im)
        χ = contribution(rel, ω, k)
        χref = contribution(refs, ω, k)
        @test_broken χ ≈ χref rtol = 1.0e-4
    end
    @test_logs (:warn, r"damped superluminal") contribution(rel, 0.9 - 0.01im, k)
    continue_lower(ys, target) = begin
        samples = [contribution(rel, 0.9 + y * im, k) for y in ys]
        V = [y^j for y in ys, j in 0:8]
        vt = [target^j for j in 0:8]
        reshape([sum((V \ [χ[i] for χ in samples]) .* vt) for i in 1:9], 3, 3)
    end
    χa = continue_lower(range(0.01, 0.15, 12), -0.05)
    χb = continue_lower(range(0.02, 0.2, 12), -0.05)
    @test_broken χa ≈ χb rtol = 1.0e-4
end
