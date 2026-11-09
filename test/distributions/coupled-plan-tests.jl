# The FixedNodeEval backend: fixed-node per-k plan for CoupledVDF on an inseparable kappa f₀.
# Oracles: the closed-form BiKappa χ (shares no quadrature with the plan; κ tails truncated by the
# box limit agreement to ~1e-6), and the exact adaptive path on the same box. Large k⊥ engages the
# product integration (Bessel weights against the coarse-column interpolant); deep damping takes the
# Landau residue far below the axis, where that interpolant — fitted on the real axis — degrades.
@testitem "CoupledVDF plan ≡ BiKappa closed form ≡ adaptive (large k⊥, k∥ ≷ 0, deep damping)" begin
    using VlasovMaxwellDispersion: plan_contribution
    using PlasmaPopulations: pdf
    bk = BiKappa(; vth_perp = 0.73, vth_para = 0.577, kappa = 3.0)
    s = prepare(NormalizedPopulation(1.0, 1.0, CoupledVDF((q, u) -> pdf(bk, q, u); para = 12.0, perp = 15.0)))
    sb = NormalizedPopulation(1.0, 1.0, bk)
    for k in (Wavevector(6.0, 0.4), Wavevector(2.0, -0.4), Wavevector(4.0, 0.0))
        pl = plan_contribution(s, k)
        @test pl isa VlasovMaxwellDispersion.TablePlan
        for ω in (1.3 + 0.1im, 1.3 + 0.0im, 1.3 - 0.3im, 0.5 - 1.5im)
            χ = pl(ω)
            @test χ ≈ contribution(sb, ω, k) rtol = 1.0e-5
            @test χ ≈ contribution(s, ω, k; rtol = 1.0e-10) rtol = 1.0e-8
        end
    end
    @test !(plan_contribution(s, Wavevector(1.0, 0.4); backend = AdaptiveEval()) isa VlasovMaxwellDispersion.TablePlan)
end

# A Landau pole GRAZING a fixed node weights it wₗ=uwₗ/(uₗ−ζ) ~ 1/δ: whatever the peeled value
# φ(ζ) disagrees with the node table by is amplified by uw/|δ|. Large k⊥ is where they can disagree
# (the table is a coarse-column interpolant there), on the real axis and just below it (crossed).
@testitem "CoupledVDF plan: Landau pole grazing a quadrature node at large k⊥" begin
    using VlasovMaxwellDispersion: plan_contribution
    s = prepare(NormalizedPopulation(-1.0, 1.0, CoupledVDF((v, u) -> exp(-(u^2 + v^2 + 0.6u * v)); para = 8.0, perp = 6.0)))
    k = Wavevector(6.0, 1.0)                     # k∥=1 ⇒ ζ_{n=0} = ω, so δ = u★ − ω exactly
    pl = plan_contribution(s, k)
    u★ = pl.para.un[findfirst(>(0.6), pl.para.un)]
    for δ in (0.0, 1.0e-16, 1.0e-12, 1.0e-8, 1.0e-6, 1.0e-4), ω in (u★ + δ + 0.0im, u★ - δ * im)
        @test pl(ω) ≈ contribution(s, ω, k; rtol = 1.0e-10) rtol = 1.0e-8
    end
end

@testitem "FixedNodeEval knobs and relativistic fallback" begin
    using VlasovMaxwellDispersion: plan_contribution
    using PlasmaPopulations: pdf
    f0(q, u) = exp(-(sqrt(q^2 + u^2) - 1.0)^2 / 0.3^2)
    d = CoupledVDF(f0; para = 2.0, perp = 2.0)
    s = NormalizedPopulation(-1.0, 1.0, d)
    k = Wavevector(0.6, 0.4)
    ω = 1.2 + 0.05im
    # a coarser grid converges to the same χ, whether built per k or once in `prepare`
    coarse = FixedNodeEval(order = 8, nprobe = 65, tol = 1.0e-7)
    χ = plan_contribution(s, k)(ω)
    @test plan_contribution(s, k; backend = coarse)(ω) ≈ χ rtol = 1.0e-6
    sc = prepare(s; backend = coarse)
    @test sc.vdf.cache.table.backend == coarse
    @test plan_contribution(sc, k; backend = coarse)(ω) ≈ χ rtol = 1.0e-6
    @test isnothing(prepare(s; backend = AdaptiveEval()).vdf.cache.table)
    # a problem prepares its populations, given as a tuple too, with its own backend
    sep = NormalizedPopulation(1.0, 1.0, SeparableVDF(Gaussian(0.1) ⊗ Kappa(0.1, 4); para = 1.0, perp = 1.0))
    prob = DispersionProblem((s, sep), (0.5 - 0.1im, 1.5 + 0.1im), k; backend = coarse)
    prepped = VlasovMaxwellDispersion.pointsolver(prob, AAA(); prob.kwargs...).prob.plasma.populations
    @test prepped[1].vdf.cache.table.backend == coarse
    @test prepped[2].vdf isa VlasovMaxwellDispersion.PreparedVDF

    mx = BiMaxwellian(; vth_perp = 0.05, vth_para = 0.1)
    rel = CoupledVDF((q, u) -> pdf(mx, q, u); para = 0.6, perp = 0.6, regime = Relativistic())
    srel = NormalizedPopulation(-1.0, 1.0, rel)
    krel = Wavevector(0.7, 0.4)
    pl = plan_contribution(srel, krel; backend = FixedNodeEval())
    @test pl(0.3 + 0.02im) ≈ contribution(srel, 0.3 + 0.02im, krel) rtol = 1.0e-6
end
