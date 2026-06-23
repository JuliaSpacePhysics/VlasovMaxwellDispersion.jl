@testitem "Maxwellian -> cold limit" begin
    # vth -> 0 must reproduce the Stix cold ε to ~vth^2 accuracy (no resonance
    # crossing: Omega large compared to omega so zeta stays far from the pole).
    Omega_e, Pi2 = -1836.0, 0.1
    k = Wavenumber(0.2, 0.3)
    omega = 0.5 + 0.0im

    epsilon_cold = dielectric(Species(Omega_e, Pi2, ColdVDF()), omega, k)
    for vth in (1.0e-2, 1.0e-3, 1.0e-4)
        epsilon_hot = dielectric(Species(Omega_e, Pi2, Maxwellian(vth)), omega, k)
        @test maximum(abs.(epsilon_hot .- epsilon_cold)) < 10 * vth^2
    end
end

@testitem "Langmuir wave: electrostatic_det matches standard Z dispersion" begin
    using VlasovMaxwellDispersion: Z, VlasovMaxwellDispersion as VM

    # 1 + (1/k^2 lambdaD^2)[1 + zeta Z(zeta)] = 0, lambdaD^2 = vth^2/(2 Pi2)
    # (vth = sqrt(2T/m) here, vs. textbook vt = sqrt(T/m) => lambdaD = vt/wp).
    # Omega_e large & purely parallel k decouples cyclotron structure so
    # electrostatic_det reduces to the unmagnetized chi_zz Langmuir branch.
    Omega_e, Pi2, vth = -1836.0, 1.0, 1.0
    k = Wavenumber(0.0, 1.0)
    kz = VM.para(k)
    lambdaD2 = vth^2 / (2 * Pi2)

    pl = Plasma(Species(Omega_e, Pi2, Maxwellian(vth)))
    f_es(omega) = electrostatic_det(pl, omega, k) / VM.abs2(k)
    f_std(omega) = 1 + (1 / (kz^2 * lambdaD2)) * (1 + (omega / (kz * vth)) * Z(omega / (kz * vth)))

    omega_BG = sqrt(Complex(Pi2 + 3 * kz^2 * vth^2 / 2))  # Bohm-Gross seed
    r_es = VM.muller(f_es, omega_BG - 1.0e-3, omega_BG + 1.0e-3, omega_BG + 1.0e-3im)
    r_std = VM.muller(f_std, omega_BG - 1.0e-3, omega_BG + 1.0e-3, omega_BG + 1.0e-3im)

    @test !isnan(r_es)
    @test isapprox(r_es, r_std; atol = 1.0e-9)        # k*lambdaD ~ 0.71: damped root, exact Z-formula match
    @test imag(r_es) < -0.3                        # genuine Landau damping, not the undamped real root
end

@testitem "Langmuir wave: Bohm-Gross + weak Landau damping" begin
    using VlasovMaxwellDispersion: VlasovMaxwellDispersion as VM

    # Moderate k*lambdaD (~0.25): Bohm-Gross real part good to ~1%, damping
    # small but resolvable (gamma/omega_r ~ -2e-3), cross-checked against the
    # same standard Z-function dispersion used above.
    Omega_e, Pi2, vth, kz = -1836.0, 1.0, 0.5, 0.7
    k = Wavenumber(0.0, kz)
    lambdaD2 = vth^2 / (2 * Pi2)
    klD = sqrt(kz^2 * lambdaD2)

    pl = Plasma(Species(Omega_e, Pi2, Maxwellian(vth)))
    f_es(omega) = electrostatic_det(pl, omega, k) / VM.abs2(k)
    f_std(omega) = 1 + (1 / (kz^2 * lambdaD2)) * (1 + (omega / (kz * vth)) * Z(omega / (kz * vth)))

    omega_BG = sqrt(Complex(Pi2 + 3 * kz^2 * vth^2 / 2))
    seed = omega_BG - 0.01im  # nudge off the real axis so Muller finds the damped root
    r_es = VM.muller(f_es, seed - 1.0e-3, seed + 1.0e-3, seed + 1.0e-3im)
    r_std = VM.muller(f_std, seed - 1.0e-3, seed + 1.0e-3, seed + 1.0e-3im)

    @test klD < 0.3
    @test isapprox(r_es, r_std; atol = 1.0e-9)
    @test isapprox(real(r_es), real(omega_BG); rtol = 0.02)  # Bohm-Gross, weak-damping regime
    @test imag(r_es) < 0  # damped (Landau)
    @test imag(r_es) / real(r_es) < -1.0e-3                  # damping resolvable, not numerical noise
end
