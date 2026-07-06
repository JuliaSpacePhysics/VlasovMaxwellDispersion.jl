# Exactly perpendicular propagation kz = 0:
# Non-relativistic parallel kernel 1/(ω − kz·u − nΩ) loses its u-pole,
# Relativistic paths keep u-poles (ωγ = nΩ) and are exercised at Im ω > 0 only.

@testitem "kz=0 Bernstein" begin
    using VlasovMaxwellDispersion: contribution, muller
    using SpecialFunctions: besselix

    vth, Pi2, kperp = 0.5, 2.0, 1.3
    s = NormalizedSpecies(1.0, Pi2, Maxwellian(vth))
    λ = (kperp * vth)^2 / 2
    # k∥→0 limit of 1 + (1/k²λD²)[1 + ζ₀ Σ Λ_n Z(ζ_n)]:  Λ_n = I_n e^{-λ}, Ω = 1
    εl(ω) = 1 + (2Pi2 / (kperp * vth)^2) * (1 - sum(besselix(abs(n), λ) * ω / (ω - n) for n in -80:80))

    @testset "Maxwellian εxx matches the I_n(λ) series" begin
        for ω in (1.37, 1.37 + 0.21im, 0.83 - 0.15im)   # damped ω needs no continuation at kz=0
            χ = contribution(s, ω, Wavenumber(kperp, 0.0))
            @test 1 + χ[1, 1] ≈ εl(complex(ω)) rtol = 1.0e-6
        end
    end

    @testset "Bernstein root via electrostatic_de" begin
        k = Wavenumber(kperp, 0.0)
        seed = sqrt(complex(1 + Pi2))   # upper-hybrid branch, between n=1 and n=2
        r = muller(ω -> electrostatic_det(s, ω, k), seed - 1.0e-3, seed, seed + 1.0e-3im)
        rref = muller(εl, seed - 1.0e-3, seed, seed + 1.0e-3im)
        @test r ≈ rref rtol = 1.0e-8
        @test abs(imag(r)) < 1.0e-10   # Bernstein modes are undamped at kz=0
    end
end

@testitem "kz=0 cross-path agreement: Maxwellian vs quadrature paths" begin
    using VlasovMaxwellDispersion: contribution, Newberger

    sp(vdf) = NormalizedSpecies(1.0, 1.0, vdf)
    f0(q, u) = exp(-(q^2 + u^2) / 0.16)
    sep = SeparableVDF(q -> exp(-q^2 / 0.16), u -> exp(-u^2 / 0.16); para = (-3.0, 3.0), perp = (0.0, 3.0))
    cpl = CoupledVDF(f0; para = (-3.0, 3.0), perp = (0.0, 3.0))
    k0 = Wavenumber(0.7, 0.0)
    for ω in (1.37, 1.4 + 0.2im, 1.4 - 0.1im)
        ref = contribution(sp(Maxwellian(0.4)), ω, k0)
        rel(x) = maximum(abs.(x .- ref)) / maximum(abs.(ref))
        @test rel(contribution(sp(sep), ω, k0)) < 1.0e-5
        @test rel(contribution(sp(cpl), ω, k0)) < 1.0e-5
        @test rel(contribution(sp(cpl), ω, k0; closure = Newberger())) < 1.0e-5
    end
    # even f∥, no drift ⇒ χxz = χyz = 0 exactly on the analytic path
    χ = contribution(sp(Maxwellian(0.4)), 1.4 + 0.1im, k0)
    scale = maximum(abs.(χ))
    @test all(abs.((χ[1, 3], χ[3, 1], χ[2, 3], χ[3, 2])) .≤ 1.0e-12 * scale)
end

@testitem "kz=0 kappa paths: analytic vs quadrature twins" begin
    using VlasovMaxwellDispersion: contribution

    sp(vdf) = NormalizedSpecies(1.0, 1.0, vdf)
    k0 = Wavenumber(0.7, 0.0)
    for κ in (2, 2.5)   # integer-residue and ₂F₁ branches share the kz=0 moment formula
        a = (κ - 1.5) * 0.4^2
        pbk = ProductBiKappa(vth_para = 0.4, kappa_para = κ)
        sep = SeparableVDF(
            v -> (1 + v^2 / a)^(-(κ + 1)), u -> (1 + u^2 / a)^(-(κ + 1));
            para = (-8.0, 8.0), perp = (0.0, 8.0)
        )
        bk = BiKappa(vth_para = 0.4, vth_perp = 0.5, kappa = κ)
        fbk(q, u) = (1 + u^2 / bk.a_para + q^2 / bk.a_perp)^(-(κ + 1))
        cplk = CoupledVDF(fbk; para = (-8.0, 8.0), perp = (0.0, 8.0))
        for ω in (1.4 + 0.2im, 1.4 - 0.1im)
            χp = contribution(sp(pbk), ω, k0)
            @test maximum(abs.(contribution(sp(sep), ω, k0) .- χp)) / maximum(abs.(χp)) < 1.0e-4
            χb = contribution(sp(bk), ω, k0)
            @test maximum(abs.(contribution(sp(cplk), ω, k0) .- χb)) / maximum(abs.(χb)) < 1.0e-4
        end
    end
end

@testitem "kz=0 GridVDF fast path ≡ coupled path on the same fit" begin
    using VlasovMaxwellDispersion: contribution

    vperp = range(0, 2.5, 40)
    vpar = range(-2.5, 2.5, 60)
    f = [exp(-(q^2 + u^2) / 0.16) for q in vperp, u in vpar]
    g = GridVDF(collect(vperp), collect(vpar), f)
    sp = NormalizedSpecies(1.0, 1.0, g)
    spc = NormalizedSpecies(1.0, 1.0, g.coupled)
    for ω in (0.5 + 0.1im, 0.5 - 0.05im)
        a = contribution(sp, ω, Wavenumber(0.2, 0.0))
        b = contribution(spc, ω, Wavenumber(0.2, 0.0))
        @test maximum(abs.(a .- b)) / maximum(abs.(a)) < 1.0e-5
    end
end

@testitem "kz→0 continuity: χ(±kz→0) → χ(0)" begin
    using VlasovMaxwellDispersion: contribution

    sp(vdf) = NormalizedSpecies(1.0, 1.0, vdf)
    vdfs = (
        Maxwellian(0.4),
        Maxwellian(vth_para = 0.4, vd = 0.3),          # drift: exercises the odd moments
        ProductBiKappa(vth_para = 0.4, kappa_para = 2.5),
        BiKappa(vth_para = 0.4, kappa = 2),   # integer κ: ₂F₁ path takes ~2min/call at ζ~10⁴
    )
    ω = 1.4 + 0.2im
    for vdf in vdfs, kz in (1.0e-4, -1.0e-4)
        χ0 = contribution(sp(vdf), ω, Wavenumber(0.7, 0.0))
        χε = contribution(sp(vdf), ω, Wavenumber(0.7, kz))
        @test χε ≈ χ0 rtol = 1.0e-4
    end
end

@testitem "kz=0 relativistic (Im ω > 0): CoupledVDF vs MaxwellJuttner Swanson" begin
    using VlasovMaxwellDispersion: contribution

    sp(vdf) = NormalizedSpecies(1.0, 1.0, vdf)
    mj(q, u) = exp(-8.0 * (sqrt(1 + q^2 + u^2) - 1))
    cvdf = CoupledVDF(mj; para = (-2.0, 2.0), perp = (0.0, 2.0), regime = Relativistic())
    for ω in (0.5 + 0.05im, 1.3 + 0.1im)   # γ-resonances inside support, off-axis at Im ω>0
        a0 = contribution(sp(cvdf), ω, Wavenumber(0.3, 0.0))
        aε = contribution(sp(cvdf), ω, Wavenumber(0.3, 1.0e-4))
        b = contribution(sp(MaxwellJuttner(mu = 8.0)), ω, Wavenumber(0.3, 0.0))
        # kz=0 is the limit of the backend's own kz→0 values …
        @test a0 ≈ aε rtol = 1.0e-3
        # … and cross-checks Swanson at the fixed-GL box path's baseline accuracy
        @test a0 ≈ b rtol = 5.0e-3
    end
end
