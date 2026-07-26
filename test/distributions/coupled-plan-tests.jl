# The FixedNodeEval backend: fixed-node per-k plan for CoupledVDF, exact analog of LowRankVDF
# without the surrogate. Must reproduce the adaptive (AdaptiveEval) path bit-for-physics.

@testitem "CoupledVDF FixedNodeEval ≡ AdaptiveEval ≡ direct" begin
    f0(q, u) = exp(-(sqrt(q^2 + u^2) - 1.0)^2 / 0.3^2)     # inseparable ring, genuinely 2-D
    d = CoupledVDF(f0; para = 2.0, perp = 2.0)
    s = NormalizedSpecies(-1.0, 1.0, d)
    # ω list includes crossed/damped poles; k∥ = 0 exercises the harmonic-independent branch
    for (k, ωs) in ((Wavenumber(0.6, 0.4), (1.2 + 0.05im, 0.7 + 0.02im, 1.4 - 0.05im, 1.4 - 0.3im)),
            (Wavenumber(0.6, 0.0), (1.2 + 0.05im, 0.7 + 0.02im, 1.5 + 0.3im)))
        plp = plan_contribution(s, k; backend = FixedNodeEval())
        pla = plan_contribution(s, k; backend = AdaptiveEval())
        @test plp isa VlasovMaxwellDispersion.CoupledPlan
        @test !(pla isa VlasovMaxwellDispersion.CoupledPlan)   # falls back to GenericKPlan
        for ω in ωs
            χ = contribution(s, ω, k)
            @test plp(ω) ≈ χ rtol = 1.0e-6
            @test plp(ω) ≈ pla(ω) rtol = 1.0e-6
        end
    end
end

# A Landau pole GRAZING a fixed node is worse than landing on it: wₗ=uwₗ/(uₗ−ζ) ~ 1/δ, so
# both Σₗwₗφₗ and φ(ζ)·Σₗwₗ reach ~uw/|δ| and must cancel back to O(1) — a loss of
# eps·uw/|δ|. Only the exact hit is removable by inspection, so the near miss needs its own
# guard. Requires Im ζ ≈ 0 too, since Im δ = −Im ζ alone keeps |δ| large.
@testitem "CoupledVDF plan: Landau pole grazing a quadrature node" begin
    f0(q, u) = exp(-(sqrt(q^2 + u^2) - 1.0)^2 / 0.3^2)
    d = CoupledVDF(f0; para = 2.0, perp = 2.0)
    s = NormalizedSpecies(-1.0, 1.0, d)
    k = Wavenumber(0.6, 1.0)                     # k∥=1 ⇒ ζ_{n=0} = ω, so δ = u★ − ω exactly
    pl = plan_contribution(s, k; backend = FixedNodeEval())
    u★ = pl.un[50]
    for ε in (0.0, 1.0e-16, 1.0e-12, 1.0e-8, 1.0e-4)
        ω = complex(u★ + ε, 0.0)
        @test all(isfinite, pl(ω))
        @test pl(ω) ≈ contribution(s, ω, k) rtol = 1.0e-6
    end
end

@testitem "FixedNodeEval is the default and threads through DispersionProblem" begin
    f0(q, u) = exp(-(sqrt(q^2 + u^2) - 1.0)^2 / 0.3^2)
    ring = CoupledVDF(f0; para = 2.0, perp = 2.0)
    plasma = (NormalizedSpecies(1.0, 0.5, ring), NormalizedSpecies(-1.0, 1.0, Maxwellian(0.2)))
    k = Wavenumber(0.6, 0.4)
    prob_p = DispersionProblem(plasma, 1.0 + 0.05im, k)                          # default backend
    prob_a = DispersionProblem(plasma, 1.0 + 0.05im, k; backend = AdaptiveEval())
    @test prob_p.backend isa FixedNodeEval
    for ω in (1.0 + 0.05im, 1.6 - 0.1im)
        @test prob_p.f(ω) ≈ prob_a.f(ω) rtol = 1.0e-6      # same det via both backends
    end
end

@testitem "FixedNodeEval knobs and relativistic fallback" begin
    f0(q, u) = exp(-(sqrt(q^2 + u^2) - 1.0)^2 / 0.3^2)
    d = CoupledVDF(f0; para = 2.0, perp = 2.0)
    s = NormalizedSpecies(-1.0, 1.0, d)
    k = Wavenumber(0.6, 0.4)
    ω = 1.2 + 0.05im
    # a coarser grid is still the same plan type and converges to the same χ
    @test plan_contribution(s, k; backend = FixedNodeEval(order = 8, nprobe = 65))(ω) ≈
        plan_contribution(s, k; backend = FixedNodeEval())(ω) rtol = 1.0e-6

    mx = Maxwellian(vth=(0.05, 0.1))
    rel = CoupledVDF(mx; para = 0.6, perp = 0.6, regime = Relativistic())
    srel = NormalizedSpecies(-1.0, 1.0, rel)
    krel = Wavenumber(0.7, 0.4)
    pl = plan_contribution(srel, krel; backend = FixedNodeEval())
    @test !(pl isa VlasovMaxwellDispersion.CoupledPlan)     # no relativistic plan yet
    @test pl(0.3 + 0.02im) ≈ contribution(srel, 0.3 + 0.02im, krel) rtol = 1.0e-6
end
