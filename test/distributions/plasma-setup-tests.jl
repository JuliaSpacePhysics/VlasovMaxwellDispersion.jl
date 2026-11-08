# Pins the physical → normalized map: (B0, n, vdf) → (Omega, Pi2, vdf in c), and `scales`.

@testitem "physical setup reproduces bo_case1's hand normalization" begin
    using Unitful
    using VlasovMaxwellDispersion: solver_vdf
    c0, qe, mp, me, eps0, mu0 = 2.99792458e8, 1.602176634e-19, 1.67262192369e-27, 9.1093837015e-31, 8.8541878128e-12, 1.25663706212e-6
    B0, n = 0.1, 5.0e19
    wcp = qe * B0 / mp
    Pi2(m) = n * qe^2 / (eps0 * m) / wcp^2
    vth(T, m) = sqrt(2qe * T / m) / c0

    p, e = Proton(), Electron()
    plasma = Plasma(
        Population(p, BiKappa; T_perp = 993.367u"eV", T_para = 1986.734u"eV", kappa = 6, n = 5.0e19u"m^-3"),
        Population(e, BiMaxwellian; T = 496.683, n = 5.0e19);
        B0 = 0.1u"T",
    )
    sp, se = NormalizedPlasma(plasma).populations
    @test sp.Omega == 1.0
    @test sp.Pi2 ≈ Pi2(mp)
    @test se.Omega ≈ -mp / me
    @test se.Pi2 ≈ Pi2(me)
    @test sp.vdf == solver_vdf(BiKappa(; vth_perp = vth(993.367, mp), vth_para = vth(1986.734, mp), kappa = 6), 1)
    @test se.vdf == solver_vdf(BiMaxwellian(; vth = vth(496.683, me)), 1)

    s = scales(plasma)
    @test s.Omega_ref ≈ wcp
    @test s.vA ≈ B0 / sqrt(mu0 * n * (mp + me)) / c0   # total mass density
    @test s.populations[1].wp ≈ sqrt(Pi2(mp))
    @test 1 / s.populations[1].d ≈ sqrt(Pi2(mp))           # case-study `kunit`s
    @test s.populations[1].rho ≈ vth(993.367, mp)          # |Ω̃| = 1
    @test s.populations[2].vth_para ≈ vth(496.683, me)
end

@testitem "reference frequency: default first particle, any particle, or a frequency" begin
    using PlasmaPopulations: E_SI, MP_SI, ME_SI
    p, e = Proton(), Electron()
    n, B0 = 5.0e6, 5.0e-9
    plasma = Plasma(Population(p, Cold(); n), Population(e, Cold(); n); B0)
    @test NormalizedPlasma(plasma).populations[1].Omega == 1.0
    se = NormalizedPlasma(plasma; ref = e).populations
    @test se[2].Omega == -1.0
    @test se[1].Omega ≈ ME_SI / MP_SI
    # normalized to ωcp, a particle the plasma does not contain
    α = Particle(:He, 2)
    @test NormalizedPlasma(Plasma(Population(α, Cold(); n); B0); ref = p).populations[1].Omega ≈ 2p.m / α.m

    Ωp = E_SI * B0 / MP_SI
    only1(pl; ref) = only(NormalizedPlasma(pl; ref).populations)
    a = only1(Plasma(Population(p, Cold(); n); B0); ref = p)
    b = only1(Plasma(Population(p, Cold(); n); B0); ref = Ωp)
    @test a.Omega ≈ b.Omega && a.Pi2 ≈ b.Pi2
    # an arbitrary frequency: Omega carries B0, Pi2 does not
    s1 = only1(Plasma(Population(p, Cold(); n); B0); ref = 2.0e3)
    s2 = only1(Plasma(Population(p, Cold(); n); B0 = 2B0); ref = 2.0e3)
    @test s2.Omega ≈ 2 * s1.Omega
    @test s1.Pi2 ≈ s2.Pi2
end

@testitem "SolverVDFs pass through; models convert or are rejected" begin
    using VlasovMaxwellDispersion: Juttner, Oriented
    p = Proton()
    plasma(vdf) = Plasma(Population(p, vdf; n = 1.0e6); B0 = 1.0e-8)
    sep = Gaussian(0.1) ⊗ Gaussian(0.1)
    @test NormalizedPlasma(plasma(sep)).populations[1].vdf === sep
    @test isnothing(scales(plasma(sep)).populations[1].vth_para)
    @test NormalizedPlasma(plasma(MaxwellJuttner(; theta = 0.1))).populations[1].vdf == Juttner(10.0)
    @test NormalizedPlasma(plasma(Cold())).populations[1].vdf isa Cold
    @test_throws ArgumentError NormalizedPlasma(plasma(Cold(; vd_para = 1.0e3)))
    @test_throws ArgumentError NormalizedPlasma(plasma(BiKappa(; vth = 1.0e4, kappa = 3, vd_para = 1.0e3)))
    @test_throws ArgumentError NormalizedPlasma(plasma(ProductBiKappa(; vth = 1.0e4, kappa_perp = 3, kappa_para = 3, loss_cone = 1)))
    @test_throws ArgumentError NormalizedPlasma(plasma(Oriented(Cold(), [1.0, 0, 0])))
end

@testitem "GridVDF converts to the SplineVDF of its data" begin
    v_perp, v_para = range(0, 0.06, 41), range(-0.06, 0.06, 61)
    F = [exp(-(v^2 + u^2) / 0.01^2) for v in v_perp, u in v_para]
    s = NormalizedPopulation(-1.0, 1.0, SplineVDF(v_perp, v_para, F))
    g = NormalizedPopulation(-1.0, 1.0, GridVDF(v_perp, v_para, F))
    ω, k = 1.2 + 0.01im, Wavevector(0.5, 0.3)
    @test contribution(g, ω, k) ≈ contribution(s, ω, k)
end

@testitem "negative density subtracts a population" begin
    using VlasovMaxwellDispersion: dielectric
    p, n, B0 = Proton(), 1.0e6, 1.0e-8
    ω, k = 0.9 - 0.01im, Wavevector(100.0, 200.0)
    for vdf in (BiMaxwellian(; vth = 1.0e5), ProductBiKappa(; vth = 1.0e5, kappa_perp = 4, kappa_para = 4), BiKappa(; vth = 1.0e5, kappa = 4), Cold())
        ref = NormalizedPlasma(Plasma(Population(p, vdf; n); B0))
        split = NormalizedPlasma(Plasma(Population(p, vdf; n = 2n), Population(p, vdf; n = -n); B0))
        @test split.populations[2].Pi2 ≈ -ref.populations[1].Pi2
        @test dielectric(split, ω, k) ≈ dielectric(ref, ω, k)
    end
end

@testitem "units follow the plasma: physical Plasma is SI" begin
    using Unitful
    using VlasovMaxwellDispersion: dielectric
    p, e = Proton(), Electron()
    plasma = Plasma(Population(p, BiMaxwellian; T = 100.0, n = 5.0e6), Population(e, BiMaxwellian; T = 50.0, n = 5.0e6); B0 = 5.0e-9)
    normalized = NormalizedPlasma(plasma)
    Ω, ku = normalized.scales.Omega_ref, normalized.scales.Omega_ref / 2.99792458e8
    @test NormalizedPlasma(Plasma(plasma.populations...; B0 = [3.0e-9, 0.0, 4.0e-9])).scales.Omega_ref ≈ Ω
    SI(k) = Wavevector(k.k_perp * ku, k.k_para * ku)

    k̃ = Wavevector(0.0, 0.5 * normalized.scales.populations[1].wp)
    ω̃ = 0.65 - 0.001im
    prob = DispersionProblem(plasma, ω̃ * Ω, SI(k̃))
    @test prob.omega0 == ω̃ * Ω && prob.k == SI(k̃)
    point = solve(prob)
    @test point.prob === prob
    @test point.omega ≈ Ω * solve(DispersionProblem(normalized, ω̃, k̃)).omega
    @test VlasovMaxwellDispersion.solve!(VlasovMaxwellDispersion.init(prob)).omega == point.omega
    @test point.resid < 1.0e-10
    @test VlasovMaxwellDispersion.residual(point.prob, point.omega) ≈ point.resid atol = 1.0e-10

    path̃ = [Wavevector(0.0, f * k̃.k_para) for f in (0.9, 1.0, 1.1)]
    tracked = solve(DispersionProblem(plasma, Seed(ω̃ * Ω, SI(path̃[2])), SI.(path̃)))
    @test tracked.omega ≈ Ω .* solve(DispersionProblem(normalized, Seed(ω̃, path̃[2]), path̃)).omega

    box̃ = (0.0 - 0.3im, 1.1 + 0.05im)
    sweep̃ = CartesianSweep(; k_para = [0.9, 1.0] .* k̃.k_para)
    sweep = CartesianSweep(; k_para = [0.9, 1.0] .* k̃.k_para .* ku)
    survey, reference = solve(DispersionProblem(plasma, box̃ .* Ω, sweep)), solve(DispersionProblem(normalized, box̃, sweep̃))
    same(a, b) = isnan.(a) == isnan.(b) && filter(isfinite, a) ≈ filter(isfinite, b)   # NaN marks missing branch points
    @test length(survey) == length(reference) > 0
    @test all(zip(survey, reference)) do (b, r)
        same(b.omega, Ω .* r.omega) && all(x -> x[1].k_para ≈ x[2].k_para * ku, zip(b.k, r.k))
    end
    # The growth margin stays relative to Ω_ref rather than becoming an absolute rate
    @test map(isgrowing, survey) == map(isgrowing, reference)
    polar = AngleSweep(; k = [0.9, 1.0] .* k̃.k_para .* ku, theta = 0.0)
    @test same(solve(DispersionProblem(plasma, box̃ .* Ω, polar))[1].omega, survey[1].omega)

    @test dielectric(plasma, ω̃ * Ω, SI(k̃)) ≈ dielectric(normalized, ω̃, k̃)
    @test_throws ArgumentError VlasovMaxwellDispersion.DispersionFunction(plasma, SI(k̃))
end
