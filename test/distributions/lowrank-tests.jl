# Rank-revealing separable surrogate of a coupled f₀. Oracles: closed forms that the
# surrogate never sees (bi-Maxwellian, BiKappa) and the exact adaptive coupled path.

@testitem "LowRankVDF: separable f₀ is rank 1 and reproduces the bi-Maxwellian" begin
    using PlasmaPopulations: pdf
    using LinearAlgebra: rank
    vthp, vthq = 0.9, 1.2
    mx = BiMaxwellian(; vth_perp = vthq, vth_para = vthp)
    lr = LowRankVDF((q, u) -> pdf(mx, q, u); para = 10vthp, perp = 10vthq)
    @test rank(lr) == 1
    for kperp in (0.1, 1.5), ω in (1.3 - 0.05im, 1.3 + 0.05im, 1.3 - 2.0im)
        k = Wavevector(kperp, 0.4)
        χl = contribution(NormalizedPopulation(-1.0, 0.5, lr), ω, k)
        χm = contribution(NormalizedPopulation(-1.0, 0.5, mx), ω, k)
        @test χl ≈ χm rtol = 1.0e-6
    end
end

# BiKappa is genuinely coupled with algebraic tails, and has an independent closed form.
@testitem "LowRankVDF: BiKappa vs its closed form, incl. damped ω" begin
    using VlasovMaxwellDispersion: plan_contribution
    using PlasmaPopulations: pdf
    using LinearAlgebra: rank
    vthz, vthp, κ = 0.2, 0.3, 3.0
    bk = BiKappa(; vth_perp = vthp, vth_para = vthz, kappa = κ)
    lr = LowRankVDF((q, u) -> pdf(bk, q, u); para = 20vthz, perp = 20vthp, rtol = 1.0e-10)
    @test 5 <= rank(lr) <= 25
    sb = NormalizedPopulation(1.0, 1.0, bk)
    sl = prepare(NormalizedPopulation(1.0, 1.0, lr))
    for (kperp, kz) in ((2.0, 1.0), (8.0, 0.1))
        k = Wavevector(kperp, kz)
        plan = plan_contribution(sl, k)
        for ω in (1.3 + 0.05im, 1.3 - 0.02im, 1.3 - 0.1im)
            @test plan(ω) ≈ contribution(sb, ω, k) rtol = 1.0e-4
        end
    end
end

# The whole point of the skeleton: its parallel factors are slices of the TRUE f₀, so the
# Landau residue is exact. A fitted surrogate (e.g. SplineVDF's spline) fails this test.
@testitem "LowRankVDF: inseparable f₀ tracks the exact coupled path under damping" begin
    g0(v, u) = exp(-(u^2 + v^2 + 0.6u * v))
    kw = (para = 8.0, perp = 6.0)
    cpl = prepare(NormalizedPopulation(-1.0, 1.0, CoupledVDF(g0; kw...)))
    lr = prepare(NormalizedPopulation(-1.0, 1.0, LowRankVDF(g0; kw..., rtol = 1.0e-10)))
    for kperp in (0.3, 1.2)
        k = Wavevector(kperp, 0.4)
        for ω in (1.2 + 0.05im, 1.2 - 0.05im, 1.2 - 0.4im)
            @test contribution(lr, ω, k) ≈ contribution(cpl, ω, k) rtol = 1.0e-5
        end
    end
end

@testitem "LowRankVDF: k∥=0 (perpendicular) path" begin
    g0(v, u) = exp(-(u^2 + v^2 + 0.6u * v))
    kw = (para = 8.0, perp = 6.0)
    cpl = NormalizedPopulation(-1.0, 1.0, CoupledVDF(g0; kw...))
    lr = prepare(NormalizedPopulation(-1.0, 1.0, LowRankVDF(g0; kw..., rtol = 1.0e-10)))
    k = Wavevector(0.8, 0.0)
    @test contribution(lr, 1.2 + 0.05im, k) ≈ contribution(cpl, 1.2 + 0.05im, k) rtol = 1.0e-5
end

# Regression for the Landau-Cauchy branches of `_lr_cauchy`. The exact adaptive CoupledVDF is
# the oracle on every path: the FAR Neumann branch's crossed-pole residue, the near-field
# `_peel` gate that avoids catastrophic cancellation when bₛ(ζ) grows off-axis, and the
# signed-zero real pole. Rank-1 Gaussian, as in the reviewer's counterexamples.
@testitem "LowRankVDF: Landau-Cauchy edge cases match the exact coupled path" begin
    g0(v, u) = exp(-(u^2 + v^2))
    kw = (para = 8.0, perp = 6.0)           # U=8 ⇒ far/near split at |ζ|=16
    cpl = prepare(NormalizedPopulation(1.0, 1.0, CoupledVDF(g0; kw...)))
    lr = prepare(NormalizedPopulation(1.0, 1.0, LowRankVDF(g0; kw..., rtol = 1.0e-10)))
    cases = [
        (Wavevector(0.4, 0.1), 1.0 - 2.0im),          # far branch, crossed n=1 pole (|ζ|=20)
        (Wavevector(0.4, 0.1), 1.0 - 2.5im),          # deeper far crossing (|ζ|=25)
        (Wavevector(0.4, 0.1), 0.6im),                # near branch, bₛ(ζ)=e³⁶ off-axis (cancellation)
        (Wavevector(0.4, 0.1), 0.4im),                # near branch, growing bₛ(ζ)
        (Wavevector(0.4, 0.5), complex(0.3, -0.0)),   # signed −0 real pole (Plemelj half-residue)
        (Wavevector(1.0, 0.2), 1.2 - 0.05im),
        (Wavevector(1.0, 0.2), 1.2 - 0.4im),          # deep damped
    ]
    for (k, ω) in cases
        @test contribution(lr, ω, k) ≈ contribution(cpl, ω, k) rtol = 1.0e-6
    end
end

# A real pole landing exactly on a fixed quadrature node is the 0*safe_inv(0)→NaN hazard;
# GRAZING one is the harder case (wₗ~1/δ ⇒ Σwₗbₛ and bₛ(ζ)Σwₗ cancel from ~uw/|δ| to O(1),
# losing eps·uw/|δ|)
@testitem "LowRankVDF: coincident node and negative-Ω refinement" begin
    using VlasovMaxwellDispersion: plan_contribution
    g0(v, u) = exp(-(u^2 + v^2))
    kw = (para = 8.0, perp = 6.0)
    cplraw, lrraw = CoupledVDF(g0; kw...), LowRankVDF(g0; kw..., rtol = 1.0e-10)

    # k∥=1 ⇒ ζ_{n=0}=ω; ω = a GL node ⇒ ζ hits that node exactly, ω+ε ⇒ grazes it.
    sc = prepare(NormalizedPopulation(1.0, 1.0, cplraw))
    sl = prepare(NormalizedPopulation(1.0, 1.0, lrraw))
    k = Wavevector(0.6, 1.0)
    plan = plan_contribution(sl, k)

    u★ = plan.pa.un[argmax(abs.(plan.pa.Bv[:, 1]))]
    for ε in (0.0, 1.0e-16, 1.0e-12, 1.0e-8, 1.0e-4)
        ω = complex(u★ + ε, 0.0)
        @test all(isfinite, plan(ω))
        @test plan(ω) ≈ contribution(sc, ω, k) rtol = 1.0e-6
    end

    # Negative Ω at k⊥=7 (a=−7): refinement must see the Bessel oscillation via |a|.
    snc = prepare(NormalizedPopulation(-1.0, 1.0, cplraw))
    snl = prepare(NormalizedPopulation(-1.0, 1.0, lrraw))
    kh, ωh = Wavevector(7.0, 0.4), 1.2 - 0.05im
    @test contribution(snl, ωh, kh) ≈ contribution(snc, ωh, kh) rtol = 1.0e-6
end

@testitem "LowRankVDF: rtol drives rank and accuracy; plan is approximate" begin
    using VlasovMaxwellDispersion: plan_contribution
    using PlasmaPopulations: pdf
    using LinearAlgebra: rank
    using VlasovMaxwellDispersion: isexact
    vthz, vthp = 0.2, 0.3
    bk = BiKappa(; vth_perp = vthp, vth_para = vthz, kappa = 3.0)
    sb = NormalizedPopulation(1.0, 1.0, bk)
    k, ω = Wavevector(4.0, 0.5), 1.3 - 0.02im
    err = map((1.0e-6, 1.0e-10)) do rtol
        lr = LowRankVDF((q, u) -> pdf(bk, q, u); para = 20vthz, perp = 20vthp, rtol)
        s = prepare(NormalizedPopulation(1.0, 1.0, lr))
        χ = contribution(s, ω, k)
        (rank(lr), maximum(abs, χ .- contribution(sb, ω, k)) / maximum(abs, contribution(sb, ω, k)))
    end
    @test err[2][1] > err[1][1]        # tighter rtol ⇒ higher rank
    @test err[2][2] < err[1][2] / 10   # ⇒ at least an order more accurate
    lr = LowRankVDF((q, u) -> pdf(bk, q, u); para = 20vthz, perp = 20vthp)
    plan = plan_contribution(prepare(NormalizedPopulation(1.0, 1.0, lr)), k)
    @test !isexact(plan)
end

@testitem "LowRankVDF: survey drops roots the surrogate cannot resolve" begin
    using PlasmaPopulations: pdf
    using VlasovMaxwellDispersion: trust_error, trusted, residual
    vthz, vthp = 0.2, 0.3
    bk = BiKappa(; vth_perp = vthp, vth_para = vthz, kappa = 3.0)
    d = LowRankVDF((q, u) -> pdf(bk, q, u); para = 20vthz, perp = 20vthp, rtol = 1.0e-10)
    # the cross is fitted on the real axis and degrades off it — sharply, in the tails first
    @test maximum(trust_error(d, u + 0.0im) for u in -4:0.5:4) < 1.0e-7
    @test maximum(trust_error(d, u - 0.3im) for u in -4:0.5:4) > 1.0e-4

    s = NormalizedPopulation(1.0, 1.0, d)
    k = Wavevector(2.0, 1.0)
    exact = NormalizedPopulation(1.0, 1.0, CoupledVDF((q, u) -> pdf(bk, q, u); para = 20vthz, perp = (0.0, 20vthp)))
    for ω in (1.3 + 0.3im, 1.3 + 0.2im)
        χl = contribution(s, ω, k)
        χe = contribution(exact, ω, k)
        @test χl ≈ χe rtol = 1.0e-6
        @test trusted(d, s, ω, k)
    end
    @testset "A damped root of the surrogate's determinant is not a root of the exact determinant" begin
        ω_fake = 1.02719 - 0.30161im
        @test !trusted(d, s, ω_fake, k)
        @test residual(s, ω_fake, k) < 1.0e-5
        @test residual(exact, ω_fake, k) > 1.0e-3
    end

    sol = solve(DispersionProblem(s, (0.2 - 0.4im, 1.5 + 0.1im), k))
    for b in sol
        @test trusted(d, s, b.omega, k)
        @test residual(exact, b.omega, k) < 1.0e-6
    end
    @test length(sol) == 2
end
