# Full magnetized EM tensor for arbitrary separable analytic f, validated by the
# Gaussian⊗Gaussian ≡ bi-Maxwellian identity (same physics, independent code path:
# generic `hilbert` parallel moments + Bessel-quadrature perp moments vs Z/Γ_n).

@testitem "SeparableVDF(Gaussian) χ matches bi-Maxwellian at oblique k" begin
    using VlasovMaxwellDispersion
    using VlasovMaxwellDispersion: contribution
    vthp, vthq = 0.9, 1.2
    sep = SeparableVDF(u -> exp(-(u / vthp)^2) / (sqrt(pi) * vthp),
                       v -> exp(-(v / vthq)^2) / (pi * vthq^2);
                       parlower=-14vthp, parupper=14vthp, perpupper=14vthq)
    mx = Maxwellian(vth_par=vthp, vth_perp=vthq)
    for (Ω, Pi2, ω, kz, kp) in ((-1.0, 0.5, 1.3 - 0.05im, 0.4, 0.3),
                                (-1.0, 0.5, 0.7 + 0.02im, 0.25, 0.6),
                                (2.0, 0.8, 1.1 - 0.1im, 0.5, 0.2))
        k = Wavenumber(kp, kz)
        χs = contribution(Species(Ω, Pi2, sep), ω, k)
        χm = contribution(Species(Ω, Pi2, mx), ω, k)
        @test maximum(abs.(χs .- χm)) / maximum(abs.(χm)) < 1e-8
    end
end

@testitem "SeparableVDF oblique dispersion root matches Maxwellian" begin
    using VlasovMaxwellDispersion
    vthp, vthq = 0.05, 0.05
    sep = SeparableVDF(u -> exp(-(u / vthp)^2) / (sqrt(pi) * vthp),
                       v -> exp(-(v / vthq)^2) / (pi * vthq^2);
                       parlower=-14vthp, parupper=14vthp, perpupper=14vthq)
    k = Wavenumber(0.2, 0.3)
    ions = Species(1.0, 1 / 1836, ColdVDF())
    ωs = solve(LocalDispersionProblem(Plasma(Species(-1.0, 1.0, sep), ions), k, 1.0 - 1e-3im)).omega
    ωm = solve(LocalDispersionProblem(Plasma(Species(-1.0, 1.0, Maxwellian(vth_par=vthp, vth_perp=vthq)), ions), k, 1.0 - 1e-3im)).omega
    @test abs(ωs - ωm) < 1e-8
end

@testitem "SeparableVDF requires oblique (kperp≠0)" begin
    using VlasovMaxwellDispersion
    using VlasovMaxwellDispersion: contribution
    sep = SeparableVDF(u -> exp(-u^2) / sqrt(pi), v -> exp(-v^2) / pi;
                       parlower=-10.0, parupper=10.0, perpupper=10.0)
    @test_throws ArgumentError contribution(Species(-1.0, 1.0, sep), 1.0 + 0im, Wavenumber(0.0, 0.5))
end

@testitem "SeparableVDF accepts a non-Gaussian f (finite χ)" begin
    using VlasovMaxwellDispersion
    using VlasovMaxwellDispersion: contribution
    # Generalized-Lorentzian (kappa-like) parallel × Gaussian perp — no closed form.
    fpar(u) = (1 + u^2 / 3)^(-2)
    sep = SeparableVDF(fpar, v -> exp(-v^2) / pi;
                       parlower=-30.0, parupper=30.0, perpupper=10.0)
    χ = contribution(Species(-1.0, 1.0, sep), 1.2 - 0.05im, Wavenumber(0.3, 0.4))
    @test all(isfinite, χ)
end
