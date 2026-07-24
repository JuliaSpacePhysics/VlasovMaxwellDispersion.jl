@testitem "Ring closures match generic SeparableVDF" begin
    ringsep(vthpar, vthperp, vd, vr) = SeparableVDF(
        Maxwellian(; vth_para = vthpar, vth_perp = vthperp, vd, vr);
        para = (vd - 8vthpar, vd + 8vthpar), perp = vr + 9 * vthperp / sqrt(2)
    )
    rbsep(vthpar, vthperp, vd, vr) = SeparableVDF(
        GaussianRing(; vth_para = vthpar, vth_perp = vthperp, vd, vr);
        para = (vd - 9vthpar, vd + 9vthpar), perp = vr + 12vthperp
    )
    cases = (
        (1.0, 0.4, 0.6, 0.1, 0.12, 0.0, 1.3 + 0.02im, 0.05),
        (-1.0, 0.3, 0.8, 0.1, 0.1, 0.05, 2.1 - 0.05im, 0.15),
        (2.0, 0.5, 0.5, 0.18, 0.18, 0.2, 2.1 - 0.05im, 0.4),
    )
    for (Ω, kz, kp, vthpar, vthperp, vd, ω, vr) in cases
        k = Wavenumber(kp, kz)
        gyro = Maxwellian(; vth_para = vthpar, vth_perp = vthperp, vd, vr)
        literal = GaussianRing(; vth_para = vthpar, vth_perp = vthperp, vd, vr)
        χ(d) = contribution(NormalizedSpecies(Ω, 1.0, d), ω, k)
        @test χ(gyro) ≈ χ(ringsep(vthpar, vthperp, vd, vr)) rtol = 1.0e-6
        @test χ(literal) ≈ χ(rbsep(vthpar, vthperp, vd, vr)) rtol = 1.0e-6
    end
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
    for (vthpar, vthperp, vd, vr) in (
        (0.1, 0.12, 0.0, 0.15),
        (0.1, 0.1, 0.05, 0.4),
    )
        k = Wavenumber(0.0, 0.4)   # k⊥=0
        χ(d) = contribution(NormalizedSpecies(Ω, 1.0, d), ω, k)
        gyro = Maxwellian(; vth_para = vthpar, vth_perp = vthperp, vd, vr)
        literal = GaussianRing(; vth_para = vthpar, vth_perp = vthperp, vd, vr)
        @test χ(gyro) ≈ χ(ringsep(vthpar, vthperp, vd, vr)) rtol = 1.0e-6
        @test χ(literal) ≈ χ(rbsep(vthpar, vthperp, vd, vr)) rtol = 1.0e-6
    end
end
