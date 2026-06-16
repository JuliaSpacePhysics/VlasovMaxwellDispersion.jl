# Most-general path: arbitrary f₀(v∥,v⊥), full magnetized EM.
# Kept fast with small k⊥ (few harmonics) — the equivalence is k-independent; 
# the general nested quadrature is the slow path by design.

@testitem "CoupledVDF(Gaussian) ≡ bi-Maxwellian (oblique)" begin
    vthp, vthq = 0.9, 1.2
    f0(u, v) = exp(-(u / vthp)^2) / (sqrt(pi) * vthp) * exp(-(v / vthq)^2) / (pi * vthq^2)
    cpl = CoupledVDF(f0; parlower=-10vthp, parupper=10vthp, perpupper=10vthq)
    mx = Maxwellian(vth_par=vthp, vth_perp=vthq)
    k = Wavenumber(0.1, 0.4)                     # small k⊥ ⇒ few harmonics ⇒ fast
    χc = contribution(Species(-1.0, 0.5, cpl), 1.3 - 0.05im, k)
    χm = contribution(Species(-1.0, 0.5, mx), 1.3 - 0.05im, k)
    @test χc ≈ χm rtol = 1e-7
end

@testitem "CoupledVDF inseparable f₀" begin
    g0(u, v) = exp(-(u^2 + v^2 + 0.6u * v))      # v∥–v⊥ coupling ⇒ not separable
    cpl = CoupledVDF(g0; parlower=-8.0, parupper=8.0, perpupper=6.0)
    χ = contribution(Species(-1.0, 1.0, cpl), 1.2 - 0.05im, Wavenumber(0.1, 0.4))
    @test all(isfinite, χ)
end

@testitem "CoupledVDF Newberger (A) ≡ HarmonicSum (B)" begin
    g0(u, v) = exp(-(u^2 + v^2 + 0.6u * v))
    cpl = CoupledVDF(g0; parlower=-8.0, parupper=8.0, perpupper=6.0)
    s = Species(-1.0, 1.0, cpl)
    ω = 1.2 + 0.05im
    for kperp in (0.3, 0.6)
        k = Wavenumber(kperp, 0.4)
        χB = contribution(s, ω, k)
        χA = contribution(s, ω, k; closure=Newberger())
        @test χA ≈ χB rtol = 1e-6
    end
end


# Regression: the relativistic harmonic path must carry the non-resonant e∥e∥
# Bernstein term 𝒳_B (derivation §5). In the weakly-relativistic limit (narrow f₀,
# γ≈1) the relativistic χ_zz must converge to the non-relativistic path, which
# folds Bernstein into m33 and is the trusted oracle.
@testitem "Relativistic CoupledVDF carries Bernstein term (anisotropic χ_zz)" begin
    ap, aq = 100.0, 400.0                        # narrow anisotropic Gaussian (p_th∥=0.1, p_th⊥=0.05)
    f0(u, w) = exp(-ap * u^2 - aq * w^2)
    dpar(u, w) = -2ap * u * f0(u, w)
    dperp(u, w) = -2aq * w * f0(u, w)
    L = 0.6
    cpl = CoupledVDF(f0; parlower=(-L), parupper=L, perpupper=L, dpar, dperp)
    ω, k = 0.3 + 0.02im, Wavenumber(0.7, 0.4)
    oracle = contribution(Species(1.0, 1.0, cpl), ω, k)[3, 3]                          # nonrel (m33 fold)
    s = Species(1.0, 1.0, cpl; regime=Relativistic())
    relB = contribution(s, ω, k)[3, 3]
    relA = contribution(s, ω, k; closure=Newberger())[3, 3]
    # Without 𝒳_B the relativistic χ_zz is off by ~3×; the residual here is the
    # genuine relativistic correction (γ−1≈5e-3), so a 2% tolerance is decisive.
    @test abs(relB - oracle) / abs(oracle) < 0.02
    @test abs(relA - oracle) / abs(oracle) < 0.02
end

@testitem "CoupledVDF requires oblique (kperp≠0)" begin
    cpl = CoupledVDF((u, v) -> exp(-(u^2 + v^2)); parlower=-8.0, parupper=8.0, perpupper=6.0)
    @test_throws ArgumentError contribution(Species(-1.0, 1.0, cpl), 1.0 + 0im, Wavenumber(0.0, 0.5))
end

@testitem "CoupledVDF Newberger (A) handles damped modes (residue extraction)" tags=[:slow] begin
    g0(u, v) = exp(-(u^2 + v^2 + 0.6u * v))
    kw = (parlower=-8.0, parupper=8.0, perpupper=6.0)
    s = Species(-1.0, 1.0, CoupledVDF(g0; kw...))
    k = Wavenumber(0.3, 0.4)
    for ω in (1.2 - 0.05im, 1.2 - 0.2im)             # damped: Im ω<0
        χB = contribution(s, ω, k)
        χA = contribution(s, ω, k; closure=Newberger())
        @test maximum(abs.(χA .- χB)) / maximum(abs.(χB)) < 1e-6
    end
end

# Relativistic evaluator A in (γ,p∥): the resonance ellipse straightens to a linear
# pole ζ_n(γ)=(ωγ−nΩ₀)/k∥, peeled by residue extraction with the Landau term — so it
# handles damped relativistic modes (Im ω<0), cross-validating B (vs Maxwell–Jüttner).
@testitem "Relativistic CoupledVDF Newberger (A) handles damped modes" begin
    μ = 40.0
    γ(u, w) = sqrt(1 + u^2 + w^2)
    f0(u, w) = exp(-μ * γ(u, w))
    dpar(u, w) = -μ * f0(u, w) * u / γ(u, w)
    dperp(u, w) = -μ * f0(u, w) * w / γ(u, w)
    L = sqrt((1 + 25 / μ)^2 - 1)
    kw = (parlower=(-L), parupper=L, perpupper=L, dpar=dpar, dperp=dperp)
    s = Species(1.0, 1.0, CoupledVDF(f0; kw...); regime=Relativistic())
    k = Wavenumber(0.7, 0.4)
    for ω in (0.3 - 0.05im, 0.3 - 0.005im)              # damped relativistic
        χB = contribution(s, ω, k)
        χA = contribution(s, ω, k; closure=Newberger())
        @test maximum(abs.(χA .- χB)) / maximum(abs.(χB)) < 1e-5
    end
end

@testitem "Relativistic CoupledVDF B finite at large k⊥ (in-range pole guard)" begin
    μ = 40.0
    γ(u, w) = sqrt(1 + u^2 + w^2)
    f0(u, w) = exp(-μ * γ(u, w))
    dpar(u, w) = -μ * f0(u, w) * u / γ(u, w)
    dperp(u, w) = -μ * f0(u, w) * w / γ(u, w)
    L = sqrt((1 + 25 / μ)^2 - 1)
    kw = (parlower=(-L), parupper=L, perpupper=L, dpar=dpar, dperp=dperp)
    s = Species(1.0, 1.0, CoupledVDF(f0; kw...); regime=Relativistic())
    ω = 0.3 - 0.05im
    for kperp in (1.2, 2.0, 3.5)                        # k⊥ρ≳1.5: off-disk poles appear
        k = Wavenumber(kperp, 0.4)
        χB = contribution(s, ω, k)                      # evaluator B (HarmonicSum)
        χA = contribution(s, ω, k; closure=Newberger()) # evaluator A (reference)
        @test all(isfinite, χB) && maximum(abs, χB) < 1e3     # no Iₙ overflow
        @test maximum(abs.(χA .- χB)) / maximum(abs.(χB)) < 1e-5
    end
end
