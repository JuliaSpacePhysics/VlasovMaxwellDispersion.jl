@testitem "BiKappa ≡ CoupledVDF (integer + non-integer κ)" begin
    ω = 1.2 - 0.05im
    for κ in (6, 4.5), kz in (0.0, 0.3)
        k = Wavenumber(0.4, kz)
        vdf = BiKappa(vth_para = 0.9, vth_perp = 1.2, kappa = κ)
        cpl = NormalizedSpecies(-1.0, 0.7, CoupledVDF(vdf; para = (-30.0, 30.0), perp = 30.0))
        bik = NormalizedSpecies(-1.0, 0.7, vdf)
        @test contribution(bik, ω, k) ≈ contribution(cpl, ω, k)
    end
end

@testitem "ProductBiKappa: anisotropic-index, → bi-Maxwellian as κ→∞" begin
    ω = 1.2 - 0.05im
    k = Wavenumber(0.4, 0.3)
    vp, vq = 0.9, 1.2
    vdf = ProductBiKappa(vth_para = vp, vth_perp = vq, kappa_para = 6, kappa_perp = 4)
    cpl = NormalizedSpecies(-1.0, 0.7, CoupledVDF(vdf; para = (-40.0, 40.0), perp = 40.0))
    pbk = NormalizedSpecies(-1.0, 0.7, vdf)
    @test contribution(pbk, ω, k) ≈ contribution(cpl, ω, k) rtol = 1.0e-6

    # κ→∞: → bi-Maxwellian (same thermal speeds), ~1/κ convergence (20→40 halves the deviation).
    χ(vdf) = contribution(NormalizedSpecies(-1.0, 0.7, vdf), ω, k)
    mx = χ(Maxwellian(vth_para = vp, vth_perp = vq))
    dev(κ) = maximum(abs.(χ(ProductBiKappa(vth_para = vp, vth_perp = vq, kappa_para = κ)) .- mx)) / maximum(abs, mx)
    @test dev(40) < 0.04
    @test dev(40) < dev(20) / 1.7

    # non-integer κ∥ (₂F₁ path) is continuous with the integer (residue) path
    pa = χ(ProductBiKappa(vth_para = vp, vth_perp = vq, kappa_para = 4, kappa_perp = 8))
    pb = χ(ProductBiKappa(vth_para = vp, vth_perp = vq, kappa_para = 4.0001, kappa_perp = 8))
    @test pa ≈ pb rtol = 1.0e-4

    @test_throws ArgumentError ProductBiKappa(vth_para = 1.0, kappa_para = 1.2, kappa_perp = 8)
end
