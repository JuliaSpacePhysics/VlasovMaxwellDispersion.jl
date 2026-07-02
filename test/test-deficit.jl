# Regression: the electron-deficit whistler — a cut VDF (Maxwellian × erfc) entered as a single
# arbitrary-VDF species and differentiated by the holomorphic-AD bridge (no hand gradient).
# Locks in (a) the AD path for erfc at the complex Landau argument, (b) the parallel whistler
# γ_max against the independent DirectR reduced-integral solver, (c) the scale-free residual at
# the root, and (d) a finite tensor at k⊥ = 0.

@testitem "deficit whistler: erfc entire-function VDF, parallel γ vs DirectR" begin
    using VlasovMaxwellDispersion: NormalizedSpecies, NormalizedPlasma
    using SpecialFunctions: erfc

    # electron-referenced units: ω in |Ω_e|, k̃ = k·R with R = ω_pe/Ω_e. Config βc=1.5, ξc=2, δ=1,
    # Δw=0.5 (paper eq-model, p=0, drift u0 ≈ 0).
    R = 100.0; MR = 1836.15267
    α = sqrt(1.5) / R; vcut = 2.0α; Δ = 0.5α; vmax = 9α
    S(vz) = 1 - 0.5 * erfc((vz + vcut) / (sqrt(2) * Δ))

    # full deficit f = f_M·S as ONE species — autodiff exercises the erfc holomorphic-AD bridge.
    e = SeparableVDF(vx -> exp(-vx^2 / α^2), vz -> exp(-vz^2 / α^2) * S(vz);
                     para = (-vmax, vmax), perp = (0.0, vmax))
    αp = sqrt(1.5 / MR) / R
    plasma = NormalizedPlasma(NormalizedSpecies(-1.0, R^2, e),
                              NormalizedSpecies(1 / MR, R^2 / MR, Maxwellian(vth_para = αp, vth_perp = αp)))

    # parallel whistler near k∥ = 0.405 d_e⁻¹ (DirectR reference: γ_max/|Ω_e| ≈ 0.0081, ω_r ≈ 0.126)
    k = Wavenumber(1e-3 * R, 0.405 * R)
    ω = solve(LocalDispersionProblem(plasma, k, 0.126 + 0.008im), Muller(maxiter = 60)).omega
    @test isfinite(ω)
    @test 0.10 < real(ω) < 0.16                           # quasi-parallel whistler band (u0≈0 shifts ω_r up slightly)
    @test 0.0072 < imag(ω) < 0.0090                       # γ_max/|Ω_e| matches DirectR (≈0.0081) to a few %
    @test dispersion_residual(plasma, ω, k) < 1e-6        # scale-free ⇒ genuine root

    # k⊥ = 0 exactly must give a finite tensor (transverse whistler limit, no degeneracy)
    @test all(isfinite, dispersion_tensor(plasma, ω, Wavenumber(0.0, 0.405 * R)))
end
