@testitem "Ring Maxwellian (vr≠0) closed form matches SeparableVDF ring" begin
    # Generic Bessel-moment quadrature of the SAME vdf
    ringsep(vthpar, vthperp, vd, vr) = SeparableVDF(
        Maxwellian(; vth_para = vthpar, vth_perp = vthperp, vd, vr);
        para = (vd - 8vthpar, vd + 8vthpar), perp = vr + 9 * vthperp / sqrt(2)
    )
    for (Ω, kz, kp, vthpar, vthperp, vd) in (
                (1.0, 0.4, 0.6, 0.1, 0.12, 0.0),
                (-1.0, 0.3, 0.8, 0.1, 0.1, 0.05),
                (2.0, 0.5, 0.5, 0.18, 0.18, 0.2),
            ),
            ω in (1.3 + 0.02im, 2.1 - 0.05im), vr in (0.05, 0.15, 0.4)
        k = Wavenumber(kp, kz)
        s1 = NormalizedSpecies(Ω, 1.0, Maxwellian(; vth_para = vthpar, vth_perp = vthperp, vd, vr))
        s2 = NormalizedSpecies(Ω, 1.0, ringsep(vthpar, vthperp, vd, vr))
        χ_fast = contribution(s1, ω, k)
        χ_quad = contribution(s2, ω, k)
        @test χ_fast ≈ χ_quad rtol = 1.0e-6
    end
    k = Wavenumber(0.6, 0.4)
    @test contribution(Maxwellian(; vth_para = 0.1, vth_perp = 0.12, vr = 0.0), 1.3 + 0.02im, k) ≈
        contribution(Maxwellian(; vth_para = 0.1, vth_perp = 0.12), 1.3 + 0.02im, k)
end

@testitem "GaussianRing (literal shifted-Gaussian, Route A) matches SeparableVDF" begin
    # literal magnitude-Gaussian perp f⊥=e^{-(v⊥-vr)²/c⊥²} × drifting f∥ (eq.(13) form).
    rbsep(vth_para, vth_perp, vdz, vdr) = SeparableVDF(
        GaussianRing(; vth_para, vth_perp, vd = vdz, vr = vdr);
        para = (vdz - 9vth_para, vdz + 9vth_para), perp = vdr + 12vth_perp
    )
    for (Ω, kz, kp, vth_para, vth_perp, vdz) in (
                (1.0, 0.4, 0.6, 0.1, 0.12, 0.0),
                (-1.0, 0.3, 0.8, 0.1, 0.1, 0.05),
                (2.0, 0.5, 0.5, 0.18, 0.18, 0.2),
            ),
            ω in (1.3 + 0.02im, 2.1 - 0.05im), vdr in (0.05, 0.15, 0.4)
        k = Wavenumber(kp, kz)
        d = GaussianRing(; vth_para, vth_perp, vd = vdz, vr = vdr)
        χA = contribution(NormalizedSpecies(Ω, 1.0, d), ω, k)
        χT = contribution(NormalizedSpecies(Ω, 1.0, rbsep(vth_para, vth_perp, vdz, vdr)), ω, k)
        @test χA ≈ χT rtol = 1.0e-6
    end
    # vr=0 reduces exactly to the bi-Maxwellian
    k = Wavenumber(0.6, 0.4)
    @test contribution(GaussianRing(; vth_para = 0.1, vth_perp = 0.12), 1.3 + 0.02im, k) ==
        contribution(Maxwellian(; vth_para = 0.1, vth_perp = 0.12), 1.3 + 0.02im, k)
end

@testitem "Parallel propagation (k⊥=0): ring energy survives β=0 fallback" begin
    # At k⊥=0 the Rₙ Bessel structure carries n/β. GaussianRing evaluates it as a genuine
    # finite moment of (J_{n-1}+J_{n+1}) (no n/β); GyroRing's Γ_n^ring closure can't, so it
    # short-circuits to the energy-matched Gaussian (⟨v⊥²⟩=vth²+vr²). Both must agree with the
    # β→0⁺ limit / quadrature ground truth.
    ringsep(vthpar, vthperp, vd, vr) = SeparableVDF(
        Maxwellian(; vth_para = vthpar, vth_perp = vthperp, vd, vr);
        para = (vd - 8vthpar, vd + 8vthpar), perp = vr + 9 * vthperp / sqrt(2)
    )
    rbsep(cpar, cperp, vdz, vdr) = SeparableVDF(
        GaussianRing(; vth_para = cpar, vth_perp = cperp, vd = vdz, vr = vdr);
        para = (vdz - 9cpar, vdz + 9cpar), perp = vdr + 12cperp
    )
    Ω, ω = 1.0, 1.3 + 0.02im
    for (vthpar, vthperp, vd) in ((0.1, 0.12, 0.0), (0.1, 0.1, 0.05)), vr in (0.15, 0.4)
        k = Wavenumber(0.0, 0.4)   # k⊥=0
        χ_gyro = contribution(NormalizedSpecies(Ω, 1.0, Maxwellian(; vth_para = vthpar, vth_perp = vthperp, vd, vr)), ω, k)
        @test χ_gyro ≈ contribution(NormalizedSpecies(Ω, 1.0, ringsep(vthpar, vthperp, vd, vr)), ω, k) rtol = 1.0e-6
        χ_lit = contribution(NormalizedSpecies(Ω, 1.0, GaussianRing(; vth_para = vthpar, vth_perp = vthperp, vd = vd, vr = vr)), ω, k)
        @test χ_lit ≈ contribution(NormalizedSpecies(Ω, 1.0, rbsep(vthpar, vthperp, vd, vr)), ω, k) rtol = 1.0e-6
    end
end
