# Most-general path: arbitrary f₀, full magnetized EM

@testitem "CoupledVDF(Gaussian) ≡ bi-Maxwellian (oblique)" begin
    vthp, vthq = 0.9, 1.2
    mx = Maxwellian(vth_para = vthp, vth_perp = vthq)
    cpl = CoupledVDF(mx; para = (-10vthp, 10vthp), perp = 10vthq)
    k = Wavenumber(0.1, 0.4)                     # small k⊥ ⇒ few harmonics ⇒ fast
    χc = contribution(NormalizedSpecies(-1.0, 0.5, cpl), 1.3 - 0.05im, k)
    χm = contribution(NormalizedSpecies(-1.0, 0.5, mx), 1.3 - 0.05im, k)
    @test χc ≈ χm
    # strongly damped ω ⇒ g(ζ)~1e13: exercises the direct/far conditioning branch per harmonic
    ωd = 1.3 - 2.0im
    χcd = contribution(NormalizedSpecies(-1.0, 0.5, cpl), ωd, k)
    χmd = contribution(NormalizedSpecies(-1.0, 0.5, mx), ωd, k)
    @test χcd ≈ χmd
end

@testitem "CoupledVDF Newberger (A) ≡ HarmonicSum (B) for inseparable f₀" begin
    g0(u, v) = exp(-(u^2 + v^2 + 0.6u * v))
    cpl = CoupledVDF(g0; para = (-8.0, 8.0), perp = 6.0)
    s = NormalizedSpecies(-1.0, 1.0, cpl)
    χ = contribution(s, 1.2 - 0.05im, Wavenumber(0.1, 0.4))
    @test all(isfinite, χ)
    ω = 1.2 + 0.05im
    for kperp in (0.0, 0.3, 0.6)
        k = Wavenumber(kperp, 0.4)
        χB = contribution(s, ω, k)
        χA = contribution(s, ω, k; closure = Newberger())
        @test χA ≈ χB
    end
end


# Regression: the relativistic harmonic path must carry the non-resonant e∥e∥
# Bernstein term 𝒳_B (derivation §5). In the weakly-relativistic limit (narrow f₀,
# γ≈1) the relativistic χ_zz must converge to the non-relativistic path, which
# folds Bernstein into m33 and is the trusted oracle.
@testitem "Relativistic CoupledVDF carries Bernstein term (anisotropic χ_zz)" begin
    mx = Maxwellian(vth_para = 0.1, vth_perp = 0.05)  # narrow anisotropic Gaussian
    L = 0.6
    cpl = CoupledVDF(mx; para = (-L, L), perp = L)
    cpl_rel = CoupledVDF(mx; para = (-L, L), perp = L, regime = Relativistic())
    ω, k = 0.3 + 0.02im, Wavenumber(0.7, 0.4)
    oracle = contribution(cpl, ω, k)[3, 3]                          # nonrel (m33 fold)
    relB = contribution(cpl_rel, ω, k)[3, 3]
    relA = contribution(cpl_rel, ω, k; closure = Newberger())[3, 3]
    # Without 𝒳_B the relativistic χ_zz is off by ~3×; the residual here is the
    # genuine relativistic correction (γ−1≈5e-3), so a 2% tolerance is decisive.
    @test abs(relB - oracle) / abs(oracle) < 0.02
    @test abs(relA - oracle) / abs(oracle) < 0.02
end

@testitem "CoupledVDF Newberger (A) handles damped modes (residue extraction)" tags = [:slow] begin
    g0(u, v) = exp(-(u^2 + v^2 + 0.6u * v))
    kw = (para = (-8.0, 8.0), perp = 6.0)
    s = NormalizedSpecies(-1.0, 1.0, CoupledVDF(g0; kw...))
    k = Wavenumber(0.3, 0.4)
    for ω in (1.2 - 0.05im, 1.2 - 0.2im)             # damped: Im ω<0
        χB = contribution(s, ω, k)
        χA = contribution(s, ω, k; closure = Newberger())
        @test χA ≈ χB rtol = 1.0e-6
    end
end


@testitem "CoupledVDF perp lower bound skips the empty ring core" begin
    vth_para, vth_perp, vr = 0.1, 0.05, 0.6       # vr/vth⊥=12 ⇒ core density ~e⁻¹⁴⁴
    f0 = Maxwellian(; vth_para, vth_perp, vr)     # callable gyro-ring density
    pk = (para = (-8vth_para, 8vth_para),)
    full = CoupledVDF(f0; pk..., perp = vr + 9vth_perp)
    cut = CoupledVDF(f0; pk..., perp = (vr - 6vth_perp, vr + 9vth_perp))
    k, ω = Wavenumber(0.6, 0.4), 1.3 + 0.02im
    @test contribution(full, ω, k) ≈ contribution(cut, ω, k) rtol = 1.0e-6
    @test contribution(full, ω, k; closure = Newberger()) ≈ contribution(cut, ω, k; closure = Newberger()) rtol = 1.0e-6
end

# --- Bounds-free (sinc/cot) path -----------------------------------------------

@testitem "CoupledVDF bounds-free ≡ bi-Maxwellian across ω regimes" begin
    mx = Maxwellian(vth_para = 0.9, vth_perp = 1.2)
    cpl = CoupledVDF(mx)                                # no bounds
    k = Wavenumber(1.0, 0.4)
    s(v) = NormalizedSpecies(-1.0, 0.5, v)
    # damped / growing / exactly-real / marginal / strongly-damped: the uniform
    # cot+i formula must hold on both sides of Im ω = 0 with no branch glitch
    for ω in (1.3 - 0.05im, 1.3 + 0.05im, 1.3 + 0.0im, 0.9 - 1.0e-10im, 1.3 - 2.0im)
        @test contribution(s(cpl), ω, k) ≈ contribution(s(mx), ω, k) rtol = 1.0e-6
    end
end

@testitem "CoupledVDF bounds-free: heavy-tail kappa (finite box truncates)" begin
    # κ=2: a (-30,30)/30 box truncates to ~2e-5 relative error; infinite bounds must
    # hit the 1e-6 target against the BiKappa closed form
    vdf = BiKappa(vth_para = 0.9, vth_perp = 1.2, kappa = 2.0)
    cpl = NormalizedSpecies(-1.0, 0.7, CoupledVDF(vdf))
    bik = NormalizedSpecies(-1.0, 0.7, vdf)
    ω, k = 1.2 - 0.05im, Wavenumber(0.4, 0.3)
    @test contribution(cpl, ω, k) ≈ contribution(bik, ω, k) rtol = 1.0e-6
end

@testitem "CoupledVDF bounds-free: drifted beam exercises the centered map" begin
    ud = 2.0
    f0(q, u) = exp(-(q^2 + (u - ud)^2))
    inf = CoupledVDF(f0)
    box = CoupledVDF(f0; para = (ud - 7.0, ud + 7.0), perp = 7.0)
    s(v) = NormalizedSpecies(-1.0, 0.5, v)
    k = Wavenumber(0.4, 0.3)
    @test inf.upar ≈ ud rtol = 1.0e-3
    for ω in (1.2 - 0.05im, 1.2 + 0.05im)
        @test contribution(s(inf), ω, k) ≈ contribution(s(box), ω, k) rtol = 1.0e-6
    end
end

@testitem "CoupledVDF bounds-free: narrow ring (GyroRing overflow-safe)" begin
    vth_para, vth_perp, vr = 0.1, 0.05, 0.6
    f0 = Maxwellian(; vth_para, vth_perp, vr)
    @test isfinite(f0(50.0, 0.0))                       # raw besseli overflowed here
    inf = CoupledVDF(f0)
    cut = CoupledVDF(f0; para = (-8vth_para, 8vth_para), perp = (vr - 6vth_perp, vr + 9vth_perp))
    k, ω = Wavenumber(0.6, 0.4), 1.3 + 0.02im
    @test contribution(inf, ω, k) ≈ contribution(cut, ω, k) rtol = 1.0e-6
end

@testitem "CoupledVDF bounds validation and closure/regime requirements" begin
    g0(q, u) = exp(-(q^2 + u^2))
    @test_throws ArgumentError CoupledVDF(g0; para = (-8.0, 8.0))            # mixed
    @test_throws ArgumentError CoupledVDF(g0; regime = Relativistic())       # rel needs box
    s = NormalizedSpecies(-1.0, 1.0, CoupledVDF(g0))
    @test_throws ArgumentError contribution(s, 1.2 - 0.05im, Wavenumber(0.3, 0.4); closure = Newberger())
end
