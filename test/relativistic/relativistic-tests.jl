# Relativistic path: an isotropic Maxwell–Jüttner f₀ fed through the
# general CoupledVDF must reproduce the closed Maxwell–Jüttner (Swanson) tensor
@testitem "Relativistic CoupledVDF reproduces Maxwell–Jüttner" begin
    μ = 40.0
    L = sqrt((1 + 25 / μ)^2 - 1)
    ref = MaxwellJuttner(mu=μ)
    rel = CoupledVDF(ref; para=(-L, L), perp=L, regime=Relativistic())

    groups = (
        ((0.3 - 0.005im, 0.3 + 0.05im), (0.0, 0.3, 0.6, 0.7), 1.0e-5, 1.0e-5),
        ((0.3 - 0.05im,), (0.0, 0.3, 0.6, 0.7), 1.0e-5, sqrt(eps())),
        ((0.3 - 0.05im,), (1.2, 2.0, 3.5), 1.0e-4, sqrt(eps())),
    )
    for (ωs, kperps, rtolA, rtolB) in groups, ω in ωs, kperp in kperps
        k = Wavenumber(kperp, 0.4)
        χA = contribution(rel, ω, k; closure=Newberger())
        χB = contribution(rel, ω, k)
        χref = contribution(ref, ω, k)
        @test χA ≈ χref rtol = rtolA
        @test χB ≈ χref rtol = rtolB
    end
end

# Strongly-relativistic (μ=2) regime with resonances INSIDE the support
@testitem "Relativistic CoupledVDF matches Swanson at μ=2 (subluminal, any damping)" begin
    μ = 2.0
    ref = MaxwellJuttner(mu=μ)
    P = 8.0
    rel = CoupledVDF(ref; para=(-P, P), perp=P, regime=Relativistic())

    pperp = range(0.0, 5.0, length=61)
    ppar = range(-5.0, 5.0, length=121)
    F = [ref(w, u) for w in pperp, u in ppar]
    grel = GridVDF(pperp, ppar, F; rtol=1.0e-4, regime=Relativistic())

    k = Wavenumber(0.7, 0.4)
    for ω in (0.3 + 0.0im, 0.3 + 0.05im, 0.3 - 0.001im, 0.3 - 0.005im, 0.3 - 0.02im, 0.3 - 0.1im, 0.3 - 0.15im)
        χB = contribution(rel, ω, k)
        χA = contribution(rel, ω, k; closure=Newberger())
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

@testitem "MaxwellJuttner parallel continuation vs corrected López" begin
    pair(vdf) = (NormalizedSpecies(1.0, 1.0, vdf), NormalizedSpecies(-1.0, 1.0, vdf))
    plasma = pair(MaxwellJuttner(mu=2.0))
    # Corrected López report §4.1 fixtures: first two are inside the resonance
    # band 0 < Re(ω)²-k∥² < 1; remaining points are beyond its upper edge.
    fixtures = (
        (2.0, 2.1 - 0.2im, -2.3261717360706777 - 11.832303790883067im),
        (2.0, 2.2135943621178655 - 0.2im, -6.25087611087304 + 24.41486369191926im),
        (2.0, 2.4 - 0.2im, 13.137323801938305 + 1.3383118276264714im),
        (6.0, 6.2 - 0.1im, -0.01702911735814406 - 2.3655114970848774im),
    )
    for (kz, ω, ref) in fixtures
        f = DispersionFunction(plasma, Wavenumber(0.0, kz); mode=:L, deflate=false)
        @test f(ω) ≈ ref rtol=1.0e-10
    end
end

# Superluminal (|Re ω| > |k∥|), damped (Im ω < 0): the straight (p⊥,p∥) box is the
# wrong sheet — the apex branch point has crossed the p⊥ path. A momentum-coords VDF
# has no analytic continuation there and warns; an energy-coords VDF routes to
# transported residue cycles (docs/src/relativistic.typ), the correct germ.
@testitem "Damped superluminal apex: momentum path wrong-sheet, energy path correct" begin
    using LinearAlgebra: norm
    μ = 2.0
    ref = MaxwellJuttner(mu=μ)
    P = sqrt((1 + 16 / μ)^2 - 1)
    en = CoupledVDF((γ, u) -> exp(-μ * γ); para=(-P, P), perp=P, coords=:energy, regime=Relativistic())
    mom = CoupledVDF(ref; para=(-P, P), perp=P, regime=Relativistic())  # momentum coords: no denergy → no cycles
    
    # momentum-coords VDF at damped-superluminal ω has no analytic continuation → warns
    @test_logs (:warn, r"damped superluminal") contribution(mom, 0.7 - 0.1im, Wavenumber(0.0, 0.5))

    # near-marginal superluminal (|Re ω| > |k∥|, moderate damping): the germ and the Landau sheet differ by O(1)
    for (kz, ω) in ((0.5, 0.7 - 0.1im), (0.5, 1.2 - 0.1im))
        k = Wavenumber(0.0, kz)
        χref = contribution(ref, ω, k)
        @test contribution(en, ω, k) ≈ χref rtol = 1.0e-4
        @test_broken contribution(en, ω, k; path=:landau) ≈ χref rtol = 0.5  # momentum: wrong sheet
    end

    # deep damping (kz=2.5, ω=2.5766−0.4536im): the two sheets nearly coincide, so BOTH
    # paths match — pinning the error to the apex/wrong-sheet, not a generic bug
    let k = Wavenumber(0.0, 2.5), ω = 2.5766 - 0.4536im
        χref = contribution(ref, ω, k)
        @test contribution(en, ω, k) ≈ χref rtol = 1.0e-4
        @test contribution(en, ω, k; path=:landau) ≈ χref rtol = 1.0e-3
    end
end

# Transported residue cycles: momentum-space continuation for damped-superluminal
# oblique k⊥. Math: docs/src/relativistic.typ.
@testitem "MaxwellJuttner superluminal continuation: transported residue cycles" begin
    using LinearAlgebra: norm
    using VlasovMaxwellDispersion: _mj_cycle_contribution
    let μ = 2.0
    mj = MaxwellJuttner(mu=μ)
    sp = NormalizedSpecies(1.0, 1.0, mj)
    rel(A, B) = norm(A - B) / norm(B)
    # subluminal damped: cycles reduce to the certified (p⊥,p∥) Landau rule
    @test _mj_cycle_contribution(mj, sp, 0.3 - 0.05im, Wavenumber(0.7, 0.4)) ≈
        contribution(sp, 0.3 - 0.05im, Wavenumber(0.7, 0.4)) rtol = 1.0e-3
    # parallel superluminal: cycles reproduce the certified Swanson continuation
    for (kz, ω) in ((0.5, 0.7 - 0.1im), (0.5, 0.7 - 0.4im), (2.5, 2.5766 - 0.4536im))
        @test _mj_cycle_contribution(mj, sp, ω, Wavenumber(0.0, kz)) ≈
            contribution(sp, ω, Wavenumber(0.0, kz)) rtol = 1.0e-4
    end
    ω, kz = 2.5766 - 0.4536im, 2.5
    # Oblique light-line seam continuity (|Δ| ∝ δ) and holomorphy.
    kp = 0.5
    seam(δ) = rel(_mj_cycle_contribution(mj, sp, ω, Wavenumber(kp, real(ω) - δ)),
        _mj_cycle_contribution(mj, sp, ω, Wavenumber(kp, real(ω) + δ)))
    @test seam(0.005) < 0.05
    h = 1.0e-4
    k = Wavenumber(kp, kz)
    dre = (_mj_cycle_contribution(mj, sp, ω + h, k) - _mj_cycle_contribution(mj, sp, ω - h, k)) / (2h)
    dim = (_mj_cycle_contribution(mj, sp, ω + im * h, k) - _mj_cycle_contribution(mj, sp, ω - im * h, k)) / (2im * h)
    @test norm(dre - dim) / norm(dre) < 1.0e-4
    # A/IC root continues beyond the old k⊥ ceiling
    pair2 = (NormalizedSpecies(1.0, 1.0, mj), NormalizedSpecies(-1.0, 1.0, mj))
    sol = solve(DispersionProblem(pair2, 2.5766 - 0.4536im, Wavenumber(0.5, 2.5)))
    @test sol.resid < 1.0e-10
    @test abs(sol.omega - (2.549295 - 0.4004106im)) < 1.0e-3
    end
end

# Generic residue cycles: any relativistic CoupledVDF with an analytic
# energy-form gradient `denergy(γ,u)` continues onto the damped-superluminal
# germ automatically (routes inside the Relativistic harmonic path)
@testitem "Relativistic CoupledVDF: residue cycles via denergy" begin
    using LinearAlgebra: norm
    using VlasovMaxwellDispersion: _MJ_CYCLE_QUAD
    μ = 2.0
    ref = MaxwellJuttner(mu=μ)
    P = sqrt((1 + 16 / μ)^2 - 1)
    rel = CoupledVDF((γ, u) -> exp(-μ * γ); para=(-P, P), perp=P,
        coords=:energy, regime=Relativistic())
    ω = 2.5766 - 0.4536im
    for (kperp, tol) in ((0.0, 1.0e-4), (0.2, 5.0e-4), (0.6, 1.0e-6))
        k = Wavenumber(kperp, 2.5)
        χ = contribution(rel, ω, k; quad=_MJ_CYCLE_QUAD)
        χref = contribution(ref, ω, k)       # Swanson (k⊥=0) / MJ cycles (k⊥≠0)
        @test norm(χ - χref) / norm(χref) < tol
    end
end

@testitem "Residue cycles preserve negative-kparallel parity" begin
    using LinearAlgebra: Diagonal
    using VlasovMaxwellDispersion: _mj_cycle_contribution, NormalizedSpecies

    mj = MaxwellJuttner(mu=2.0)
    sp = NormalizedSpecies(1.0, 1.0, mj)
    ω = 2.5766 - 0.4536im
    S = Diagonal([1, 1, -1])
    χp = _mj_cycle_contribution(mj, sp, ω, Wavenumber(0.5, 2.5))
    χm = _mj_cycle_contribution(mj, sp, ω, Wavenumber(0.5, -2.5))
    @test χm ≈ S * χp * S
end

@testitem "Residue cycle q2 is stable at its light-line endpoint" begin
    using VlasovMaxwellDispersion: _cycle_endpoint, _cycle_q2

    ω, kz, N = 2.5 - 1.0e-8im, 2.5, -1.0
    γ0 = _cycle_endpoint(N, ω, kz)
    u0 = (γ0 * ω - N) / kz
    @test iszero(_cycle_q2(γ0, γ0, u0, ω, kz))
    @test isfinite(_cycle_q2(γ0 + 1.0e-8, γ0, u0, ω, kz))

    ω, N = 2.5 - 1.0e-14im, 0.0
    D = sqrt(complex(N^2 + kz^2 - ω^2))
    γ0 = _cycle_endpoint(N, ω, kz)
    u0 = (γ0 * ω - N) / kz
    δ = 1.0e-6cis(0.3)
    expected = δ * (2D / abs(kz) + (1 - (ω / kz)^2) * δ)
    @test _cycle_q2(γ0 + δ, γ0, u0, ω, kz; D, δ) ≈ expected rtol = 1.0e-12
end

@testitem "Residue cycles do not infer a zero contour from its endpoint" begin
    using LinearAlgebra: norm
    using VlasovMaxwellDispersion: _cycle_endpoint, _cycle_harmonic

    ω, Ω, kz = 0.7 - 0.1im, 1.0, 0.5
    γ0 = _cycle_endpoint(0.0, ω, kz)
    scaledUzero = (γ, u, σ) -> (γ - γ0) * exp(-γ + σ)
    @test norm(_cycle_harmonic(0, scaledUzero, 1.0, ω, Ω, kz, 0.0)) > 1.0e-8

    scaledUcov = (γ, u, σ) -> exp(-γ + σ)
    χ = _cycle_harmonic(0, scaledUcov, 1.0, ω, Ω, kz, 0.0)
    @test all(isfinite, χ)
    @test scaledUcov(800, 0, 800) == 1
end
