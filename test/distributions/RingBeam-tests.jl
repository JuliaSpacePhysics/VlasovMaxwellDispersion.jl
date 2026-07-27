@testitem "Ring closures match generic SeparableVDF" begin
    ringsep(vth, vd, vr) = SeparableVDF(
        Maxwellian(; vth, vd, vr);
        para=(vd - 8vth[2], vd + 8vth[2]), perp=vr + 9 * vth[1] / sqrt(2)
    )
    rbsep(vth, vd, vr) = SeparableVDF(
        GaussianRing(; vth, vd, vr);
        para=(vd - 9vth[2], vd + 9vth[2]), perp=vr + 12vth[1]
    )
    cases = (
        (1.0, 0.4, 0.6, (0.12, 0.1), 0.0, 1.3 + 0.02im, 0.05),
        (-1.0, 0.3, 0.8, (0.1, 0.1), 0.05, 2.1 - 0.05im, 0.15),
        (2.0, 0.5, 0.5, (0.18, 0.18), 0.2, 2.1 - 0.05im, 0.4),
    )
    for (Ω, kz, kp, vth, vd, ω, vr) in cases
        k = Wavenumber(kp, kz)
        gyro = Maxwellian(; vth, vd, vr)
        literal = GaussianRing(; vth, vd, vr)
        χ(d) = contribution(NormalizedSpecies(Ω, 1.0, d), ω, k)
        @test χ(gyro) ≈ χ(ringsep(vth, vd, vr))
        @test χ(literal) ≈ χ(rbsep(vth, vd, vr))
    end

    @testset "Parallel propagation (k⊥=0): ring energy survives β=0 fallback" begin
        # At k⊥=0 the Rₙ Bessel structure carries n/β. GaussianRing evaluates it as a genuine
        # finite moment of (J_{n-1}+J_{n+1}) (no n/β); GyroRing's Γ_n^ring closure can't, so it
        # short-circuits to the energy-matched Gaussian (⟨v⊥²⟩=vth²+vr²).
        Ω, ω = 1.0, 1.3 + 0.02im
        for (vthpar, vthperp, vd) in ((0.1, 0.12, 0.0), (0.1, 0.1, 0.05)), vr in (0.15, 0.4)
            k = Wavenumber(0.0, 0.4)   # k⊥=0
            vth = (vthperp, vthpar)
            χ(d) = contribution(NormalizedSpecies(Ω, 1.0, d), ω, k)
            gyro = Maxwellian(; vth, vd, vr)
            literal = GaussianRing(; vth, vd, vr)
            @test χ(gyro) ≈ χ(ringsep(vth, vd, vr))
            @test χ(literal) ≈ χ(rbsep(vth, vd, vr))
        end
    end
end
