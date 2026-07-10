@testitem "BiKappa ≡ CoupledVDF (integer + non-integer κ)" begin
    ω = 1.2 - 0.05im

    function test_contribution(κ, kz)
        k = Wavenumber(0.4, kz)
        vdf = BiKappa(vth_para = 0.9, vth_perp = 1.2, kappa = κ)
        cpl = NormalizedSpecies(-1.0, 0.7, CoupledVDF(vdf; para = (-30.0, 30.0), perp = 30.0))
        bik = NormalizedSpecies(-1.0, 0.7, vdf)
        contribution(bik, ω, k) ≈ contribution(cpl, ω, k)
    end

    for κ in (100, 6, 4.5), kz in (0.0, 0.3)
        @test test_contribution(κ, kz)
    end
    @test_broken test_contribution(1000, 0.0)
end

@testitem "ProductBiKappa: anisotropic-index → bi-Maxwellian as κ→∞" begin
    ω = 1.2 - 0.05im
    k = Wavenumber(0.4, 0.3)
    vp, vq = 0.9, 1.2
    vdf = ProductBiKappa(vth_para = vp, vth_perp = vq, kappa_para = 6, kappa_perp = 4)
    cpl = NormalizedSpecies(-1.0, 0.7, CoupledVDF(vdf; para = (-40.0, 40.0), perp = 40.0))
    pbk = NormalizedSpecies(-1.0, 0.7, vdf)
    @test contribution(pbk, ω, k) ≈ contribution(cpl, ω, k) rtol = 1.0e-6

    # κ→∞: → bi-Maxwellian (same thermal speeds), ~1/κ convergence
    χ(vdf) = contribution(NormalizedSpecies(-1.0, 0.7, vdf), ω, k)
    mx = χ(Maxwellian(vth_para = vp, vth_perp = vq))
    test_contribution(κ) = maximum(abs.(χ(ProductBiKappa(vth_para = vp, vth_perp = vq, kappa_para = κ)) .- mx)) / maximum(abs, mx)
    @test test_contribution(1000) ≈ 0 atol = 1.0e-3
    @test_broken test_contribution(1000.001) ≈ dev(1000)

    # non-integer κ∥ (₂F₁ path) is continuous with the integer (residue) path
    pa = χ(ProductBiKappa(vth_para = vp, vth_perp = vq, kappa_para = 4, kappa_perp = 8))
    pb = χ(ProductBiKappa(vth_para = vp, vth_perp = vq, kappa_para = 4.0001, kappa_perp = 8))
    @test pa ≈ pb rtol = 1.0e-4

    # κ∥ > 1/2 (finite T∥) and κ⊥ > 1 (finite T⊥)
    @test contribution(NormalizedSpecies(-1.0, 0.7, ProductBiKappa(vth_para = 1.0, kappa_para = 1, kappa_perp = 8)), ω, k) isa AbstractMatrix
end

@testitem "kappa: NaN ζ terminates (no reflection recursion)" begin
    # A root finder probing a wild trial ω can overflow ζ = Δ/kz into NaN
    # non-integer-κ ₂F₁ reflection branch must not recurse on abs(NaN).
    using VlasovMaxwellDispersion: _kappa_H0
    H = _kappa_H0(complex(NaN, NaN), 1.0, 5.5)
    @test isnan(real(H)) && isnan(imag(H))
    @test isnan(real(_kappa_H0(complex(Inf, NaN), 1.0, 5.5, -1)))
end
