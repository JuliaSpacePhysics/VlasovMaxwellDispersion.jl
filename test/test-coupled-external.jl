# External-reference regression for genuinely-inseparable coupled path.
# Numbers baked from one-off cross-checks:
#   - ALPS (grid+poly-fit): bi-kappa dispersion roots; our det is ~0 there.

@testitem "CoupledVDF shell ε matches LMV (Newberger) reference" tags=[:slow] begin
    using VlasovMaxwellDispersion
    # matched dimensionless setup (Ωref=Ω_proton); shell f₀=exp(-(√(u²+w²)-v_s)²/2σ²)
    k̃, θ = 0.8, 60 / 180 * π
    k = Wavenumber(k̃ * sin(θ), k̃ * cos(θ))
    ω = 0.6 - 0.006im
    vthi, vthe, vshell = 0.0014600019156727816, 0.0625625037436143, 0.004380005747018345
    σ2 = vthi^2
    f0(u, w) = exp(-(sqrt(u^2 + w^2) - vshell)^2 / (2σ2))
    dpar(u, w) = (s = sqrt(u^2 + w^2); s == 0 ? 0.0 * u : f0(u, w) * (-(s - vshell) / σ2) * (u / s))
    dprp(u, w) = (s = sqrt(u^2 + w^2); s == 0 ? 0.0 * w : f0(u, w) * (-(s - vshell) / σ2) * (w / s))
    L = vshell + 9vthi
    el = Species(-1836.2059501591832, 346868.41251994757, Maxwellian(vth_par=vthe, vth_perp=vthe))
    ish = Species(1.0, 188.90496051920377,
                  CoupledVDF(f0; parlower=-L, parupper=L, perpupper=L, dpar=dpar, dperp=dprp))
    ε = dielectric(Plasma(el, ish), ω, k)
    # LMV reference ε (row-major i,j)
    εref = ComplexF64[
        296.2155861486968-3.319922531867824im   3.7626621891157885+177.04618736508462im  -0.0025053128619864605+5.637511809770261e-5im
        -3.7626621891157885-177.04618736508462im 296.20733289587076-3.3200008059935646im  -0.01416517941366486+0.4757377172149103im
        -0.0025053128619864605+5.637511809770261e-5im 0.01416517941366486-0.4757377172149103im -966280.6019165423-19378.285104224095im]
    @test maximum(abs.(ε .- εref)) / maximum(abs.(εref)) < 1e-8
end

# Relativistic (γ,p∥) coupled path: an isotropic Maxwell–Jüttner f₀ fed through
# the general CoupledVDF must reproduce the closed Maxwell–Jüttner (Swanson)
# tensor — itself ALPS-validated. ω<Ω keeps the Swanson time-integral stable.
@testitem "Relativistic CoupledVDF reproduces Maxwell–Jüttner" begin
    using VlasovMaxwellDispersion
    kz, kp, ω, Ω, μ = 0.4, 0.7, 0.3 - 0.005im, 1.0, 40.0
    γ(u, w) = sqrt(1 + u^2 + w^2)
    f0(u, w) = exp(-μ * γ(u, w))
    dpar(u, w) = -μ * f0(u, w) * u / γ(u, w)
    dperp(u, w) = -μ * f0(u, w) * w / γ(u, w)
    L = sqrt((1 + 25 / μ)^2 - 1)                 # momentum cutoff for γmax≈1+25/μ
    rel = Species(Ω, 1.0, CoupledVDF(f0; parlower=-L, parupper=L, perpupper=L, dpar, dperp);
                  regime=Relativistic())
    ref = Species(Ω, 1.0, MaxwellJuttner(mu=μ))
    k = Wavenumber(kp, kz)
    χ = contribution(rel, ω, k)
    χref = contribution(ref, ω, k)
    @test maximum(abs.(χ .- χref)) / maximum(abs.(χref)) < 1e-5
end

# Evaluator A (Newberger closed-orbit, complex-order Bessel) in (p∥,p⊥) vs the
# closed Maxwell–Jüttner — relativistic, at Im ω>0 (A has no Landau contour). A's
# coordinates differ from evaluator B's (γ,p∥), so this also cross-checks the
# relativistic coordinate transform. A actually beats B's fixed-order GL here.
@testitem "Relativistic CoupledVDF Newberger (A) reproduces Maxwell–Jüttner" begin
    using VlasovMaxwellDispersion
    using VlasovMaxwellDispersion: contribution
    kz, kp, ω, Ω, μ = 0.4, 0.7, 0.3 + 0.05im, 1.0, 40.0     # Im ω>0 for evaluator A
    γ(u, w) = sqrt(1 + u^2 + w^2)
    f0(u, w) = exp(-μ * γ(u, w))
    dpar(u, w) = -μ * f0(u, w) * u / γ(u, w)
    dperp(u, w) = -μ * f0(u, w) * w / γ(u, w)
    L = sqrt((1 + 25 / μ)^2 - 1)
    rel = Species(Ω, 1.0, CoupledVDF(f0; parlower=-L, parupper=L, perpupper=L, dpar, dperp);
                  regime=Relativistic())
    ref = Species(Ω, 1.0, MaxwellJuttner(mu=μ))
    k = Wavenumber(kp, kz)
    @test maximum(abs.(contribution(rel, ω, k; closure=Newberger()) .- contribution(ref, ω, k))) / maximum(abs.(contribution(ref, ω, k))) < 1e-7
end

@testitem "CoupledVDF bi-kappa is a root at ALPS's dispersion root" tags=[:slow] begin
    using VlasovMaxwellDispersion
    using LinearAlgebra
    vA, me, κ = 1e-4, 5.44662e-4, 6.0
    a2 = (2κ - 3) / (2κ)
    function bikappa_species(Ω, Pi2, vth)
        s = κ * a2 * vth^2
        f(u, v) = (1 + (u^2 + v^2) / s)^(-1 - κ)
        du(u, v) = (-1 - κ) * (1 + (u^2 + v^2) / s)^(-2 - κ) * (2u / s)
        dv(u, v) = (-1 - κ) * (1 + (u^2 + v^2) / s)^(-2 - κ) * (2v / s)
        L = 12vth
        Species(Ω, Pi2, CoupledVDF(f; parlower=-L, parupper=L, perpupper=L, dpar=du, dperp=dv))
    end
    plasma = Plasma(bikappa_species(1.0, 1 / vA^2, vA),
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
    @test abs(real(ω) - real(ωref)) / abs(real(ωref)) < 2e-3
    @test abs(imag(ω) - imag(ωref)) / abs(imag(ωref)) < 5e-2
end
