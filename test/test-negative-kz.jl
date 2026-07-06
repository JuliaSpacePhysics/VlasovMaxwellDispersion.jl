# Negative k∥ (contour orientation σ = sign k∥). Parity identity for gyrotropic f₀:
#   χ(ω, k⊥, −k∥; f) = S·χ(ω, k⊥, k∥; f∘(v∥→−v∥))·S,  S = diag(1,1,−1),
# exact for every path (incl. relativistic). Damped ω exercises the σ·2πi Landau
# residue; growing ω exercises the causal-side branch choice (which was also wrong
# for k∥<0: Im ω>0 puts ζ in the lower half-plane).

@testitem "hilbert σ orientation: reflection + Schwarz identities" begin
    using VlasovMaxwellDispersion: hilbert

    g = v -> exp(-v^2) * (1 + v / 3)
    gr = t -> g(-t)
    L, U = -3.0, 3.0
    for ζ in (0.4 + 0.3im, 0.4 - 0.25im, -0.7 + 0.02im, 0.9 - 0.01im)
        # substitution v→−v maps orientation σ=−1 to the standard σ=+1 integral
        @test hilbert(g, ζ, L, U; σ = -1) ≈ -hilbert(gr, -ζ, -U, -L) rtol = 1.0e-8
        # Schwarz reflection (real-analytic g): H₋(ζ) = conj(H₊(conj ζ))
        @test hilbert(g, ζ, L, U; σ = -1) ≈ conj(hilbert(g, conj(ζ), L, U)) rtol = 1.0e-8
    end
end

@testitem "hilbert_landau_pwpoly σ orientation: Schwarz identity" begin
    using VlasovMaxwellDispersion: hilbert_landau_pwpoly

    nodes = [-2.0, -0.5, 0.3, 1.1, 2.0]
    coeffs = [[1.0, 0.2], [0.8, -0.1, 0.05], [1.2, 0.3], [0.5, -0.2, 0.1]]
    for ζ in (0.4 + 0.3im, 0.4 - 0.25im, -0.7 + 0.02im)
        @test hilbert_landau_pwpoly(coeffs, nodes, ζ, -1) ≈
            conj(hilbert_landau_pwpoly(coeffs, nodes, conj(ζ), 1)) rtol = 1.0e-12
    end
end

@testitem "χ parity under k∥ → −k∥: analytic and quadrature paths" begin
    using VlasovMaxwellDispersion: contribution, Newberger
    using LinearAlgebra

    S = Diagonal([1, 1, -1])
    flipz(X) = S * X * S
    sp(vdf) = NormalizedSpecies(1.0, 1.0, vdf)
    parity(vdf, ω; kperp = 0.2, kz = 0.6, kw...) = begin
        a = contribution(sp(vdf), ω, Wavenumber(kperp, kz); kw...)
        b = contribution(sp(vdf), ω, Wavenumber(kperp, -kz); kw...)
        maximum(abs.(b .- flipz(a))) / maximum(abs.(a))
    end

    f0(q, u) = exp(-(q^2 + u^2) / 0.16)
    sep = SeparableVDF(q -> exp(-q^2 / 0.16), u -> exp(-u^2 / 0.16); para = (-3.0, 3.0), perp = (0.0, 3.0))
    cpl = CoupledVDF(f0; para = (-3.0, 3.0), perp = (0.0, 3.0))
    for ω in (0.5 + 0.1im, 0.5 - 0.05im)   # growing AND damped (Landau residue active)
        @test parity(Maxwellian(0.4), ω) < 1.0e-12
        @test parity(ProductBiKappa(vth_para = 0.4, kappa_para = 2), ω) < 1.0e-12    # integer-M residue path
        @test parity(ProductBiKappa(vth_para = 0.4, kappa_para = 2.5), ω) < 1.0e-12  # ₂F₁ branch path
        @test parity(BiKappa(vth_para = 0.4, kappa = 2.5), ω) < 1.0e-12
        @test parity(sep, ω) < 1.0e-6
        @test parity(cpl, ω) < 1.0e-6
        @test parity(cpl, ω; closure = Newberger()) < 1.0e-6
        @test parity(ReducedVDF(u -> exp(-u^2 / 0.16); para = (-3.0, 3.0)), ω; kperp = 0.0) < 1.0e-10
    end

    # drifting f∥ is odd-asymmetric: χ(−k∥; vd) = S·χ(k∥; −vd)·S
    for ω in (0.5 + 0.1im, 0.5 - 0.05im)
        a = contribution(sp(Maxwellian(vth_para = 0.4, vd = -0.3)), ω, Wavenumber(0.2, 0.6))
        b = contribution(sp(Maxwellian(vth_para = 0.4, vd = 0.3)), ω, Wavenumber(0.2, -0.6))
        @test maximum(abs.(b .- flipz(a))) / maximum(abs.(a)) < 1.0e-12
    end
end

@testitem "χ parity under k∥ → −k∥: relativistic paths" begin
    using VlasovMaxwellDispersion: contribution, Newberger, Relativistic
    using LinearAlgebra

    S = Diagonal([1, 1, -1])
    flipz(X) = S * X * S
    sp(vdf) = NormalizedSpecies(1.0, 1.0, vdf)
    parity(vdf, ω; kperp = 0.2, kz = 0.6, kw...) = begin
        a = contribution(sp(vdf), ω, Wavenumber(kperp, kz); kw...)
        b = contribution(sp(vdf), ω, Wavenumber(kperp, -kz); kw...)
        maximum(abs.(b .- flipz(a))) / maximum(abs.(a))
    end

    mj(q, u) = exp(-8.0 * (sqrt(1 + q^2 + u^2) - 1))
    cvdf = CoupledVDF(mj; para = (-2.0, 2.0), perp = (0.0, 2.0), regime = Relativistic())
    for ω in (0.5 + 0.05im, 0.5 - 0.02im)   # subluminal (|Re ω| < k∥): damped side supported
        @test parity(cvdf, ω) < 1.0e-8
        @test parity(MaxwellJuttner(mu = 8.0), ω) < 1.0e-10
    end
    # (γ,p∥) Newberger backend supports Im ω ≥ 0 only
    @test parity(cvdf, 0.5 + 0.05im; closure = Newberger()) < 1.0e-6
end

@testitem "GridVDF at k∥<0: fast grid path ≡ coupled path on the same fit" begin
    # The NNLS fit is not exactly even in v∥, so parity holds only to fit tolerance;
    # instead validate the grid fast path's σ handling against the (independently
    # parity-tested) CoupledVDF quadrature on the identical fit.
    using VlasovMaxwellDispersion: contribution

    vperp = range(0, 2.5, 40)
    vpar = range(-2.5, 2.5, 60)
    f = [exp(-(q^2 + u^2) / 0.16) for q in vperp, u in vpar]
    g = GridVDF(collect(vperp), collect(vpar), f)
    sp = NormalizedSpecies(1.0, 1.0, g)
    spc = NormalizedSpecies(1.0, 1.0, g.coupled)
    for ω in (0.5 + 0.1im, 0.5 - 0.05im), kz in (0.6, -0.6)
        a = contribution(sp, ω, Wavenumber(0.2, kz))
        b = contribution(spc, ω, Wavenumber(0.2, kz))
        @test maximum(abs.(a .- b)) / maximum(abs.(a)) < 1.0e-5
    end
end

@testitem "Langmuir root at k∥<0 mirrors the k∥>0 root" begin
    using VlasovMaxwellDispersion: VlasovMaxwellDispersion as VM

    Omega_e, Pi2, vth, kz = -1836.0, 1.0, 0.5, 0.7
    pl = NormalizedSpecies(Omega_e, Pi2, Maxwellian(vth))
    root(kzv) = begin
        k = Wavenumber(0.0, kzv)
        f = omega -> electrostatic_det(pl, omega, k) / VM.abs2(k)
        seed = sqrt(Complex(Pi2 + 3 * kzv^2 * vth^2 / 2)) - 0.01im
        VM.muller(f, seed - 1.0e-3, seed + 1.0e-3, seed + 1.0e-3im)
    end
    rp, rm = root(kz), root(-kz)
    @test !isnan(rp) && imag(rp) < 0
    @test rm ≈ rp rtol = 1.0e-10   # even f∥ ⇒ spectrum even in k∥
end
