# Full magnetized EM tensor for arbitrary separable analytic f
# i.e.: Gaussian⊗Gaussian ≡ bi-Maxwellian identity

@testitem "SeparableVDF(Gaussian) χ matches bi-Maxwellian" begin
    vthp, vthq = 0.9, 1.2
    sep = SeparableVDF(
        v -> exp(-(v / vthq)^2) / (pi * vthq^2),
        u -> exp(-(u / vthp)^2) / (sqrt(pi) * vthp);
        parlower = -14vthp, parupper = 14vthp, perpupper = 14vthq
    )
    mx = Maxwellian(vth_par = vthp, vth_perp = vthq)
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
    sep = SeparableVDF(
        v -> exp(-(v / vthq)^2) / (pi * vthq^2),
        u -> exp(-(u / vthp)^2) / (sqrt(pi) * vthp);
        parlower = -14vthp, parupper = 14vthp, perpupper = 14vthq
    )
    k = Wavenumber(0.2, 0.3)
    ions = NormalizedSpecies(1.0, 1 / 1836, ColdVDF())
    ωs = solve(LocalDispersionProblem((NormalizedSpecies(-1.0, 1.0, sep), ions), k, 1.0 - 1.0e-3im)).omega
    ωm = solve(LocalDispersionProblem((NormalizedSpecies(-1.0, 1.0, Maxwellian(vth_par = vthp, vth_perp = vthq)), ions), k, 1.0 - 1.0e-3im)).omega
    @test ωs ≈ ωm
end

@testitem "SeparableVDF supports parallel propagation" begin
    vthp, vthq = 0.9, 1.2
    sep = SeparableVDF(
        v -> exp(-(v / vthq)^2) / (pi * vthq^2),
        u -> exp(-(u / vthp)^2) / (sqrt(pi) * vthp);
        parlower = -14vthp, parupper = 14vthp, perpupper = 14vthq
    )
    mx = Maxwellian(vth_par = vthp, vth_perp = vthq)
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

@testitem "Ring Maxwellian (vr≠0) closed form matches SeparableVDF ring" begin
    using SpecialFunctions: besseli
    # gyrotropic ring f⊥ ∝ e^{-(v²+vr²)/2σ²} I₀(vr v/σ²), σ²=vthperp²/2 — the f⊥
    # that the Maxwellian(vr=…) cold-ring⊛Maxwellian convolution closes in closed form.
    ringsep(vthpar, vthperp, vd, vr) = let σ2 = vthperp^2 / 2
        SeparableVDF(
            v -> exp(-(v^2 + vr^2) / (2σ2)) * besseli(0, vr * v / σ2),
            u -> exp(-(u - vd)^2 / vthpar^2);
            parlower = vd - 8vthpar, parupper = vd + 8vthpar, perpupper = vr + 9 * vthperp / sqrt(2)
        )
    end
    for (Ω, kz, kp, vthpar, vthperp, vd) in (
                (1.0, 0.4, 0.6, 0.1, 0.12, 0.0),
                (-1.0, 0.3, 0.8, 0.1, 0.1, 0.05),
                (2.0, 0.5, 0.5, 0.18, 0.18, 0.2),
            ),
            ω in (1.3 + 0.02im, 2.1 - 0.05im), vr in (0.05, 0.15, 0.4)
        k = Wavenumber(kp, kz)
        s1 = NormalizedSpecies(Ω, 1.0, Maxwellian(; vth_par = vthpar, vth_perp = vthperp, vd, vr))
        s2 = NormalizedSpecies(Ω, 1.0, ringsep(vthpar, vthperp, vd, vr))
        χ_fast = contribution(s1, ω, k)
        χ_quad = contribution(s2, ω, k)
        @test χ_fast ≈ χ_quad rtol = 1.0e-6
    end
    # vr=0 must be bit-identical to the plain bi-Maxwellian fast path
    k = Wavenumber(0.6, 0.4)
    @test contribution(Maxwellian(; vth_par = 0.1, vth_perp = 0.12, vr = 0.0), 1.3 + 0.02im, k) ==
        contribution(Maxwellian(; vth_par = 0.1, vth_perp = 0.12), 1.3 + 0.02im, k)
end

@testitem "RingBeam (literal shifted-Gaussian, Route A) matches SeparableVDF" begin
    # literal magnitude-Gaussian perp f⊥=e^{-(v⊥-vr)²/c⊥²} × drifting f∥ (eq.(13) form).
    rbsep(cpar, cperp, vdz, vdr) = SeparableVDF(
        v -> exp(-(v - vdr)^2 / cperp^2), u -> exp(-(u - vdz)^2 / cpar^2);
        parlower = vdz - 9cpar, parupper = vdz + 9cpar, perpupper = vdr + 12cperp
    )
    for (Ω, kz, kp, cpar, cperp, vdz) in (
                (1.0, 0.4, 0.6, 0.1, 0.12, 0.0),
                (-1.0, 0.3, 0.8, 0.1, 0.1, 0.05),
                (2.0, 0.5, 0.5, 0.18, 0.18, 0.2),
            ),
            ω in (1.3 + 0.02im, 2.1 - 0.05im), vdr in (0.05, 0.15, 0.4)
        k = Wavenumber(kp, kz)
        d = RingBeam(; vth_par = cpar, vth_perp = cperp, vd = vdz, vr = vdr)
        χA = contribution(NormalizedSpecies(Ω, 1.0, d), ω, k)
        χT = contribution(NormalizedSpecies(Ω, 1.0, rbsep(cpar, cperp, vdz, vdr)), ω, k)
        @test χA ≈ χT rtol = 1.0e-6
    end
    # vr=0 reduces exactly to the bi-Maxwellian
    k = Wavenumber(0.6, 0.4)
    @test contribution(NormalizedSpecies(1.0, 1.0, RingBeam(; vth_par = 0.1, vth_perp = 0.12)), 1.3 + 0.02im, k) ==
        contribution(NormalizedSpecies(1.0, 1.0, Maxwellian(; vth_par = 0.1, vth_perp = 0.12)), 1.3 + 0.02im, k)
end
