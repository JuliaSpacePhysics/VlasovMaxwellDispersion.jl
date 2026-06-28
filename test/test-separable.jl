# Full magnetized EM tensor for arbitrary separable analytic f
# i.e.: Gaussian⊗Gaussian ≡ bi-Maxwellian identity

@testitem "SeparableVDF(Gaussian) χ matches bi-Maxwellian" begin
    vthp, vthq = 0.9, 1.2
    mx = Maxwellian(vth_par = vthp, vth_perp = vthq)
    sep = SeparableVDF(mx; parlower = -14vthp, parupper = 14vthp, perpupper = 14vthq)
    for (Ω, Pi2, ω, kz, kp) in (
            (-1.0, 0.5, 1.3 - 0.05im, 0.4, 0.3),
            (-1.0, 0.5, 0.7 + 0.02im, 0.25, 0.6),
            (2.0, 0.8, 1.1 - 0.1im, 0.5, 0.2),
        )
        k = Wavenumber(kp, kz)
        χs = contribution(NormalizedSpecies(Ω, Pi2, sep), ω, k)
        χm = contribution(NormalizedSpecies(Ω, Pi2, mx), ω, k)
        @test χs ≈ χm rtol = 1.0e-8
    end
end

@testitem "SeparableVDF oblique dispersion root matches Maxwellian" begin
    vthp, vthq = 0.05, 0.05
    mx = Maxwellian(vth_par = vthp, vth_perp = vthq)
    sep = SeparableVDF(mx; parlower = -14vthp, parupper = 14vthp, perpupper = 14vthq)
    k = Wavenumber(0.2, 0.3)
    ions = NormalizedSpecies(1.0, 1 / 1836, ColdVDF())
    ωs = solve(LocalDispersionProblem((NormalizedSpecies(-1.0, 1.0, sep), ions), k, 1.0 - 1.0e-3im)).omega
    ωm = solve(LocalDispersionProblem((NormalizedSpecies(-1.0, 1.0, mx), ions), k, 1.0 - 1.0e-3im)).omega
    @test ωs ≈ ωm
end

@testitem "SeparableVDF supports parallel propagation" begin
    vthp, vthq = 0.9, 1.2
    mx = Maxwellian(vth_par = vthp, vth_perp = vthq)
    sep = SeparableVDF(mx; parlower = -14vthp, parupper = 14vthp, perpupper = 14vthq)
    k = Wavenumber(0.0, 0.4)
    χs = contribution(NormalizedSpecies(-1.0, 0.5, sep), 1.3 - 0.05im, k)
    χm = contribution(NormalizedSpecies(-1.0, 0.5, mx), 1.3 - 0.05im, k)
    @test all(isfinite, χs)
    @test χs ≈ χm
    @test abs(χs[1, 3]) < 1.0e-12 && abs(χs[2, 3]) < 1.0e-12  # transverse/parallel decouple
end

@testitem "SeparableVDF accepts a non-Gaussian f (finite χ)" begin
    # Generalized-Lorentzian (kappa-like) parallel × Gaussian perp — no closed form.
    fpar(u) = (1 + u^2 / 3)^(-2)
    sep = SeparableVDF(
        v -> exp(-v^2) / pi, fpar;
        parlower = -30.0, parupper = 30.0, perpupper = 10.0
    )
    χ = contribution(NormalizedSpecies(-1.0, 1.0, sep), 1.2 - 0.05im, Wavenumber(0.3, 0.4))
    @test all(isfinite, χ)
end
