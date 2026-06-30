# Most-general path: arbitrary f₀(v∥,v⊥), full magnetized EM.
# Kept fast with small k⊥ (few harmonics) — the equivalence is k-independent;
# the general nested quadrature is the slow path by design.

@testitem "CoupledVDF(Gaussian) ≡ bi-Maxwellian (oblique)" begin
    vthp, vthq = 0.9, 1.2
    mx = Maxwellian(vth_para = vthp, vth_perp = vthq)
    cpl = CoupledVDF(mx; para = (-10vthp, 10vthp), perp = 10vthq)
    k = Wavenumber(0.1, 0.4)                     # small k⊥ ⇒ few harmonics ⇒ fast
    χc = contribution(NormalizedSpecies(-1.0, 0.5, cpl), 1.3 - 0.05im, k)
    χm = contribution(NormalizedSpecies(-1.0, 0.5, mx), 1.3 - 0.05im, k)
    @test χc ≈ χm rtol = 1.0e-7
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
        @test χA ≈ χB rtol = 1.0e-6
    end
end

# Relativistic (γ,p∥) coupled path: an isotropic Maxwell–Jüttner f₀ fed through
# the general CoupledVDF must reproduce the closed Maxwell–Jüttner (Swanson)
# tensor — itself ALPS-validated. ω<Ω keeps the Swanson time-integral stable.
@testitem "Relativistic CoupledVDF reproduces Maxwell–Jüttner" begin
    μ = 40.0
    L = sqrt((1 + 25 / μ)^2 - 1)
    ref = MaxwellJuttner(mu = μ)
    rel = CoupledVDF(ref; para = (-L, L), perp = L, regime = Relativistic())

    for ω in (0.3 - 0.005im, 0.3 + 0.05im), kperp in (0.0, 0.3, 0.6)
        k = Wavenumber(kperp, 0.4)
        χA = contribution(rel, ω, k; closure = Newberger())
        χB = contribution(rel, ω, k)
        χref = contribution(ref, ω, k)
        @test χA ≈ χref rtol = 1.0e-5
        @test χB ≈ χref rtol = 1.0e-5
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

# Relativistic evaluator A in (γ,p∥): the resonance ellipse straightens to a linear
# pole ζ_n(γ)=(ωγ−nΩ₀)/k∥, peeled by residue extraction with the Landau term — so it
# handles damped relativistic modes (Im ω<0), cross-validating B (vs Maxwell–Jüttner).
@testitem "Relativistic CoupledVDF Newberger (A) handles damped modes" begin
    μ = 40.0
    γ(w, u) = sqrt(1 + u^2 + w^2)
    f0(w, u) = exp(-μ * γ(w, u))
    L = sqrt((1 + 25 / μ)^2 - 1)
    kw = (para = (-L, L), perp = L)
    vdf = CoupledVDF(f0; kw..., regime = Relativistic())
    k = Wavenumber(0.7, 0.4)
    for ω in (0.3 - 0.05im, 0.3 - 0.005im)              # damped relativistic
        χB = contribution(vdf, ω, k)
        χA = contribution(vdf, ω, k; closure = Newberger())
        @test maximum(abs.(χA .- χB)) / maximum(abs.(χB)) < 1.0e-5
    end
end

@testitem "Relativistic CoupledVDF B finite at large k⊥ (in-range pole guard)" begin
    μ = 40.0
    γ(w, u) = sqrt(1 + u^2 + w^2)
    f0(w, u) = exp(-μ * γ(w, u))
    L = sqrt((1 + 25 / μ)^2 - 1)
    kw = (para = (-L, L), perp = L)
    s = CoupledVDF(f0; kw..., regime = Relativistic())
    ω = 0.3 - 0.05im
    for kperp in (1.2, 2.0, 3.5)                        # k⊥ρ≳1.5: off-disk poles appear
        k = Wavenumber(kperp, 0.4)
        χB = contribution(s, ω, k)                      # evaluator B (HarmonicSum)
        χA = contribution(s, ω, k; closure = Newberger()) # evaluator A (reference)
        @test all(isfinite, χB) && maximum(abs, χB) < 1.0e3     # no Iₙ overflow
        @test maximum(abs.(χA .- χB)) / maximum(abs.(χB)) < 1.0e-5
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
