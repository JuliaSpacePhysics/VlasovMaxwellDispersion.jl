@testitem "BiKappa ≡ CoupledVDF (integer + non-integer κ)" begin
    using PlasmaPopulations: pdf
    function test_contribution(κ, kz)
        ω = 1.2 - 0.05im
        k = Wavevector(0.4, kz)
        vdf = BiKappa(; vth_perp = 1.2, vth_para = 0.9, kappa = κ)
        cpl = NormalizedPopulation(-1.0, 0.7, CoupledVDF((q, u) -> pdf(vdf, q, u); para=30.0, perp=30.0))
        bik = NormalizedPopulation(-1.0, 0.7, vdf)
        contribution(bik, ω, k) ≈ contribution(cpl, ω, k)
    end

    for κ in (100, 6, 4.5), kz in (0.0, 0.3)
        @test test_contribution(κ, kz)
    end
    @test_broken test_contribution(1000, 0.0)
end

@testitem "ProductBiKappa matches CoupledVDF and Maxwellian limit" begin
    using PlasmaPopulations: pdf
    ω = 1.2 - 0.05im
    k = Wavevector(0.4, 0.3)
    vp, vq = 0.9, 1.2
    for kpar in (6, 4.5)
        vdf = ProductBiKappa(; vth_perp = vq, vth_para = vp, kappa_perp = 4, kappa_para = kpar)
        cpl = NormalizedPopulation(-1.0, 0.7, CoupledVDF((q, u) -> pdf(vdf, q, u); para=40.0, perp=40.0))
        pbk = NormalizedPopulation(-1.0, 0.7, vdf)
        @test contribution(pbk, ω, k) ≈ contribution(cpl, ω, k) rtol = 1.0e-6
    end

    # κ→∞: → bi-Maxwellian (same thermal speeds), ~1/κ convergence
    χ(vdf) = contribution(NormalizedPopulation(-1.0, 0.7, vdf), ω, k)
    mx = χ(BiMaxwellian(; vth_perp = vq, vth_para = vp))
    test_contribution(κ) = maximum(abs.(χ(ProductBiKappa(; vth_perp = vq, vth_para = vp, kappa_perp = κ, kappa_para = κ)) .- mx)) / maximum(abs, mx)
    @test test_contribution(1000) ≈ 0 atol = 1.0e-3

    # non-integer κ∥ (₂F₁ path) is continuous with the integer (residue) path
    pa = χ(ProductBiKappa(; vth_perp = vq, vth_para = vp, kappa_perp = 8, kappa_para = 4))
    pb = χ(ProductBiKappa(; vth_perp = vq, vth_para = vp, kappa_perp = 8, kappa_para = 4.0001))
    @test pa ≈ pb rtol = 1.0e-4
end

@testitem "kappa: NaN ζ terminates (no reflection recursion)" begin
    # A root finder probing a wild trial ω can overflow ζ = Δ/kz into NaN
    # non-integer-κ ₂F₁ reflection branch must not recurse on abs(NaN).
    using VlasovMaxwellDispersion: _kappa_H0
    H = _kappa_H0(complex(NaN, NaN), 1.0, 5.5)
    @test isnan(real(H)) && isnan(imag(H))
    @test isnan(real(_kappa_H0(complex(Inf, NaN), 1.0, 5.5, -1)))
end
