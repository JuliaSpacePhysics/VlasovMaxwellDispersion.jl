# Pins the (particle, n, B0) → (Omega, Pi2) map.
# Physical types (Particle/Species/Plasma) come from PlasmaBase.

@testitem "Species→NormalizedSpecies" begin
    V = VlasovMaxwellDispersion.PlasmaBase
    e, p = Electron(), Proton()
    B0 = 5.0e-9
    # electron relative to proton: q ratio −1, mass ratio m_p/m_e ≈ 1836
    se = NormalizedSpecies(Species(e, ColdVDF(); n = 1.0), B0, p)
    @test se.Omega ≈ -1836.0 rtol = 1.0e-3
    # self-referenced proton: Omega = 1 exactly; Omega is B-free (same at any B0)
    sp = NormalizedSpecies(Species(p, ColdVDF(); n = 1.0), B0, p)
    @test sp.Omega == 1.0
    @test NormalizedSpecies(Species(p, ColdVDF(); n = 1.0), 3B0).Omega == 1.0
    # alpha: Z=2, m=4m_p ⇒ Omega = 2/4 = 1/2 vs proton; independent of B
    sa = NormalizedSpecies(Species(Particle(; z = 2, A = 4), ColdVDF(); n = 1.0), B0, p)
    @test sa.Omega ≈ 0.5
end


@testitem "accessor interface on physical types" begin
    using VlasovMaxwellDispersion.PlasmaBase
    C = PlasmaBase
    s = Species(Proton(), ColdVDF(); n = 5.0e6)
    @test charge(s) == C.E_SI
    @test mass(s) == C.MP_SI
    @test number_density(s) == 5.0e6
    @test distribution(s) isa ColdVDF
    @test particle(s) === Proton() || charge(particle(s)) == C.E_SI
    pl = Plasma(s; B0 = 5.0e-9)
    @test magnetic_field(pl) == 5.0e-9
    @test length(species(pl)) == 1
end

@testitem "NormalizedPlasma fixes the reference: default first, overridable" begin
    using VlasovMaxwellDispersion: NormalizedPlasma
    using VlasovMaxwellDispersion.PlasmaBase
    C = PlasmaBase
    p, e = Proton(), Electron()
    B0 = 5.0e-9
    phys = Plasma(Species(p, ColdVDF(); n = 5.0e6), Species(e, ColdVDF(); n = 5.0e6); B0)
    # default ref = first species' particle (proton)
    np = NormalizedPlasma(phys)
    @test first(np.species).Omega == 1.0
    # explicit ref = electron: electron self-ref ⇒ −1
    npe = NormalizedPlasma(phys; ref = e)
    @test last(npe.species).Omega == -1.0
    @test first(npe.species).Omega ≈ C.ME_SI / C.MP_SI
end

@testitem "Ω_ref need not be a gyrofrequency: frequency ref carries B0 in Omega" begin
    using VlasovMaxwellDispersion.PlasmaBase
    C = PlasmaBase
    p = Proton()
    n, B0 = 5.0e6, 5.0e-9
    Ωp = C.E_SI * B0 / C.MP_SI            # proton gyrofrequency
    # frequency ref equal to the proton gyrofreq must match particle-ref proton
    a = NormalizedSpecies(Species(p, ColdVDF(); n), B0)
    b = NormalizedSpecies(Species(p, ColdVDF(); n), B0, Ωp)
    @test a.Omega ≈ b.Omega
    @test a.Pi2 ≈ b.Pi2
    # normalize to an ARBITRARY frequency (e.g. ω_pe): Omega is no longer B-free
    Ωref = 2.0e3
    s1 = NormalizedSpecies(Species(p, ColdVDF(); n), B0, Ωref)
    s2 = NormalizedSpecies(Species(p, ColdVDF(); n), 2B0, Ωref)
    @test s1.Omega ≈ C.E_SI * B0 / C.MP_SI / Ωref
    @test s2.Omega ≈ 2 * s1.Omega          # carries B0
    @test s1.Pi2 ≈ s2.Pi2                   # Pi2 = (ω_ps/Ω_ref)², B0-free at fixed Ω_ref
end

# # TODO dielectric for ω_phy and k_phy
# @testitem "physically-built plasma reproduces a hand-normalized one" begin
#     using VlasovMaxwellDispersion
#     using VlasovMaxwellDispersion.PlasmaBase
#     C = PlasmaBase
#     n, B0 = 5.0e6, 5.0e-9
#     p, e = Proton(), Electron()
#     phys = Plasma(Species(p, ColdVDF(); n), Species(e, ColdVDF(); n); B0))
#     # hand map: Ω_e = −m_p/m_e, Pi2_e = Pi2_p·(m_p/m_e) (same n, ω_pe²/ω_pp² = m_p/m_e)
#     Pi2p = (sqrt(n * C.MP_SI / C.EPS0_SI) / B0)^2
#     mr = C.MP_SI / C.ME_SI
#     hand = NormalizedSpecies(1.0, Pi2p, ColdVDF()), NormalizedSpecies(-mr, Pi2p * mr, ColdVDF())
#     ω, k = 1.3 - 0.05im, Wavenumber(0.2, 0.5)
#     @test dielectric(phys, ω, k) ≈ dielectric(hand, ω, k) rtol = 1e-10
# end
