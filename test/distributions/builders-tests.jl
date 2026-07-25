# Pins the (particle, n, B0) → (Omega, Pi2) map.
# Physical types (Particle/Species/Plasma) come from PlasmaBase.

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

@testitem "Plasma fixes the reference: default first, overridable" begin
    using VlasovMaxwellDispersion: NormalizedPlasma
    using VlasovMaxwellDispersion.PlasmaBase
    C = PlasmaBase
    p, e = Proton(), Electron()
    B0 = 5.0e-9
    ss = (Species(p, ColdVDF(); n = 5.0e6), Species(e, ColdVDF(); n = 5.0e6))
    # default ref = first species' particle (proton)
    np = NormalizedPlasma(Plasma(ss...; B0))
    @test first(np.species).Omega == 1.0
    # explicit ref = electron: electron self-ref ⇒ −1
    npe = NormalizedPlasma(Plasma(ss...; B0, ref = e))
    @test last(npe.species).Omega == -1.0
    @test first(npe.species).Omega ≈ C.ME_SI / C.MP_SI
end

@testitem "Ω_ref need not be a gyrofrequency: frequency ref carries B0 in Omega" begin
    using VlasovMaxwellDispersion.PlasmaBase
    C = PlasmaBase
    p = Proton()
    n, B0 = 5.0e6, 5.0e-9
    Ωp = C.E_SI * B0 / C.MP_SI            # proton gyrofrequency
    only1(pl) = only(NormalizedPlasma(pl).species)
    # frequency ref equal to the proton gyrofreq must match particle-ref proton
    a = only1(Plasma(Species(p, ColdVDF(); n); B0))
    b = only1(Plasma(Species(p, ColdVDF(); n); B0, ref = Ωp))
    @test a.Omega ≈ b.Omega
    @test a.Pi2 ≈ b.Pi2
    # normalize to an ARBITRARY frequency (e.g. ω_pe): Omega is no longer B-free
    Ωref = 2.0e3
    s1 = only1(Plasma(Species(p, ColdVDF(); n); B0, ref = Ωref))
    s2 = only1(Plasma(Species(p, ColdVDF(); n); B0 = 2B0, ref = Ωref))
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
