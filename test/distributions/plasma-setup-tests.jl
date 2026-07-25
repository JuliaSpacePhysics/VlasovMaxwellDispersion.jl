# Pins the physical setup layer: (B0, n, T) → (Omega, Pi2, vth), the `scales` a paper's
# k-unit needs, and VDF spec resolution.

@testitem "physical setup reproduces bo_case1's hand normalization" begin
    using Unitful
    c0, qe, mp, me, eps0, mu0 = 2.99792458e8, 1.602176634e-19, 1.67262192369e-27, 9.1093837015e-31, 8.8541878128e-12, 1.25663706212e-6
    B0, n = 0.1, 5.0e19
    wcp = qe * B0 / mp
    Pi2(m) = n * qe^2 / (eps0 * m) / wcp^2
    vth(T, m) = sqrt(2qe * T / m) / c0

    pl = NormalizedPlasma(
        Plasma(
            Species(Proton(), sc -> BiKappa(sc; kappa=6); n=5.0e19u"m^-3", T=(1986.734, 993.367) .* u"eV"),
            Species(Electron(); n=5.0e19u"m^-3", T=496.683u"eV"),
            B0=0.1u"T",
        )
    )
    p, e = pl.species
    @test p.Omega == 1.0
    @test p.Pi2 ≈ Pi2(mp)
    @test e.Omega ≈ -mp / me
    @test e.Pi2 ≈ Pi2(me)
    @test p.vdf == BiKappa(vth_para=vth(1986.734, mp), vth_perp=vth(993.367, mp), kappa=6)
    @test e.vdf == Maxwellian(vth(496.683, me))

    s = scales(pl)
    @test s.Omega_ref ≈ wcp
    @test s.vA ≈ B0 / sqrt(mu0 * n * (mp + me)) / c0   # total mass density
    @test s.species[1].wp ≈ sqrt(Pi2(mp))
    @test 1 / s.species[1].d ≈ sqrt(Pi2(mp))           # case-study `kunit`s
    @test s.species[1].rho ≈ vth(993.367, mp)          # |Ω̃| = 1

    # bare numbers are eV / m⁻³ / T — same plasma, no Unitful
    bare = NormalizedPlasma(Plasma(Species(Proton(); n=5.0e19, T=1986.734), B0=0.1))
    @test bare.species[1].Pi2 ≈ Pi2(mp)
    @test bare.species[1].vdf == Maxwellian(vth(1986.734, mp))
end

@testitem "thermal inputs, VDF resolution, and a ref outside the plasma" begin
    using Unitful
    c0, qe, mp, eps0, mu0 = 2.99792458e8, 1.602176634e-19, 1.67262192369e-27, 8.8541878128e-12, 1.25663706212e-6

    # β∥ = 1 ⇒ vth∥ = B/√(μ₀ n mₛ), the species' OWN Alfvén speed
    plb = NormalizedPlasma(Plasma(Species(Proton(); n=1.0e6, beta=1.0), Species(Electron(); n=1.0e6, beta=1.0), B0=5.0e-9))
    @test scales(plb).species[1].vth ≈ 5.0e-9 / sqrt(mu0 * 1.0e6 * mp) / c0
    # vth in c or Unitful; MaxwellJuttner reads mu = 2/vth²
    @test NormalizedPlasma(Plasma(Species(Electron(), MaxwellJuttner; n=1.0e6, vth=sqrt(2/10)), B0=5.0e-9)).species[1].vdf.mu ≈ 10.0
    @test scales(NormalizedPlasma(Plasma(Species(Electron(); n=1.0e6, vth=0.1c0*u"m/s"), B0=5.0e-9))).species[1].vth ≈ 0.1

    # a callable spec sees its own thermal speeds alongside the plasma-wide vA
    pls = NormalizedPlasma(
        Plasma(Species(Proton(), sc -> Maxwellian(sc.vth + sc.vA); n=1.0e6, T=10), Species(Electron(); n=1.0e6, T=10), B0=1.0e-8)
    )
    ss = scales(pls)
    @test pls.species[1].vdf == Maxwellian(ss.species[1].vth + ss.vA)

    # ice_alpha: a D–e–α plasma normalized to ωcp, a particle it never contains, with the
    # paper's k-unit 1/λp built from a density no single species carries
    pli = NormalizedPlasma(
        Plasma(
            Species(Particle(z=2, A=4), sc -> GaussianRing(sc; vr=0.045); n=1.0e16, T=1000),
            Species(Particle(z=1, A=2); n=9.98e18, T=1000),
            Species(Electron(); n=1.0e19, T=1000),
            B0=2.1, ref=Proton(),
        )
    )
    @test pli.species[1].Omega ≈ 0.5   # Ω_α = Ω_D = ωcp/2
    @test plasma_gyro_ratio(1.0e19, mass(Proton()), 2.1) ≈ sqrt(1.0e19 * qe^2 / (eps0 * mp)) / (qe * 2.1 / mp)

    # no VDF and no thermal input ⇒ cold fluid; a spec with no thermal input ⇒ error
    @test NormalizedPlasma(Plasma(Species(Proton(); n=1.0e6), B0=1.0e-8)).species[1].vdf isa ColdVDF
    @test_throws ArgumentError Species(Proton(); n=1.0e6, T=10, beta=1.0)
end
