# External-reference regression for genuinely-inseparable coupled path.
# Numbers baked from one-off cross-checks:
#   - ALPS (grid+poly-fit): bi-kappa dispersion roots; our det is ~0 there.

@testitem "CoupledVDF shell ε matches LMV (Newberger) reference" begin
    using VlasovMaxwellDispersion
    # matched dimensionless setup (Ωref=Ω_proton); shell f₀=exp(-(√(u²+w²)-v_s)²/2σ²)
    k̃, θ = 0.8, 60 / 180 * π
    k = Wavenumber(k̃ * sin(θ), k̃ * cos(θ))
    ω = 0.6 - 0.006im
    vthi, vthe, vshell = 0.0014600019156727816, 0.0625625037436143, 0.004380005747018345
    σ2 = vthi^2
    f0(w, u) = exp(-(sqrt(u^2 + w^2) - vshell)^2 / (2σ2))
    L = vshell + 9vthi
    el = NormalizedSpecies(-1836.2059501591832, 346868.41251994757, Maxwellian(vth_par=vthe, vth_perp=vthe))
    ish = NormalizedSpecies(1.0, 188.90496051920377,
        CoupledVDF(f0; parlower=(-L), parupper=L, perpupper=L))
    ε = dielectric((el, ish), ω, k)
    # LMV reference ε (row-major i,j)
    εref = ComplexF64[
        296.2155861486968-3.319922531867824im 3.7626621891157885+177.04618736508462im -0.0025053128619864605+5.637511809770261e-5im
        -3.7626621891157885-177.04618736508462im 296.20733289587076-3.3200008059935646im -0.01416517941366486+0.4757377172149103im
        -0.0025053128619864605+5.637511809770261e-5im 0.01416517941366486-0.4757377172149103im -966280.6019165423-19378.285104224095im]
    @test ε ≈ εref rtol = 1e-8
end

@testitem "CoupledVDF bi-kappa is a root at ALPS's dispersion root" begin
    using LinearAlgebra
    vA, me, κ = 1e-4, 5.44662e-4, 6.0
    a2 = (2κ - 3) / (2κ)
    function bikappa_species(Ω, Pi2, vth)
        s = κ * a2 * vth^2
        f(v, u) = (1 + (u^2 + v^2) / s)^(-1 - κ)
        L = 12vth
        NormalizedSpecies(Ω, Pi2, CoupledVDF(f; parlower=(-L), parupper=L, perpupper=L))
    end
    plasma = (bikappa_species(1.0, 1 / vA^2, vA),
        bikappa_species(-1 / me, 1 / (me * vA^2), vA / sqrt(me)))
    # ALPS bi-kappa (κ=6) root_1, first scan point
    k = Wavenumber(1e-3 / vA, 0.03 / vA)         # k̃ = k_ALPS/vA = (300, 10)
    ωref = 0.029311 - 9.9693e-6im
    # Solve from the ALPS seed instead of only checking |det| there: the prior
    # |det|<2e-3 spot-check was fragile (this mode is near-marginal, so the scaled
    # det at ωref is already ~1e-3 — it barely separated a true root from noise and
    # never pinned the location). Polishing to the actual root is far stronger: it
    # lands a genuine zero (ndet→1e-14) AND recovers BOTH ALPS numbers — frequency
    # to 1e-3 and, unusually for a near-marginal mode, the damping rate to ~1%.
    ω = solve(LocalDispersionProblem(plasma, k, ωref)).omega
    D = dispersion_tensor(plasma, ω, k)
    ndet = abs(det(D)) / prod(norm(D[i, :]) for i in 1:3)  # scale-invariant residual
    @test ndet < 1e-9                                       # genuine isolated root
    @test real(ω) ≈ real(ωref) rtol = 2e-3
    @test imag(ω) ≈ imag(ωref) rtol = 5e-2
end
