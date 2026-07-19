@testitem "TensorReduction: selectors, validation, Longitudinal" begin
    using VlasovMaxwellDispersion: TensorReduction, check_reduction, FullDet, Circular, Longitudinal, Ordinary, Extraordinary
    @test_throws ArgumentError TensorReduction(:Z)
    @test_throws ArgumentError Circular(2)
    cold = (NormalizedSpecies(1.0, 1.0, ColdVDF()), NormalizedSpecies(-1.0, 1836.0, ColdVDF()))
    # eager geometry validation at problem construction, not mid-sweep
    prob = DispersionProblem(cold, 0.5 + 0im, [Wavenumber(0.0, kz) for kz in 0.3:0.1:0.5]; mode=:L)
    @test prob.mode === Circular(+1)
    @test_throws ArgumentError DispersionFunction(cold, Wavenumber(0.1, 0.5); mode=:L)
    # rounded axis points (k⊥ = k sinθ at θ ≈ 0 or π) still qualify as parallel
    @test DispersionProblem(cold, 0.5 + 0im, AngleSweep(1.0, [Float64(π)]); mode=:L).mode === Circular(+1)
    # Longitudinal is valid at any k: exact P factor at k⊥=0, electrostatic
    # approximation obliquely (same zeros as electrostatic_det)
    ko = Wavenumber(0.4, 0.5)
    fLong = DispersionFunction(cold, ko; mode=Longitudinal(), deflate=false)
    ω = 0.7 + 0.02im
    @test abs2(ko) * fLong(ω) ≈ electrostatic_det(cold, ω, ko) rtol = 1.0e-12
end

@testitem "Parallel factorization: det = R·L·P, convention, degeneracy" begin
    cold = (NormalizedSpecies(1.0, 1.0, ColdVDF()), NormalizedSpecies(-1.0, 1836.0, ColdVDF()))
    k0 = Wavenumber(0.0, 0.5)
    fL, fR, fP, fD = map(
        m -> DispersionFunction(cold, k0; mode=m),
        (:L, :R, :P, :det)
    )
    ω = 0.3 + 0.01im
    @test fD(ω) ≈ fR(ω) * fL(ω) * fP(ω)
    # L resonates with the positive species (ω → +Ω): |L| ≫ |R| near ω = 1
    @test abs(fL(0.999)) > 100 * abs(fR(0.999))
    # equal-mass pair plasma: R and L exactly degenerate ⇒ every det root is
    # a double zero; the factors have simple zeros
    pair2 = (NormalizedSpecies(1.0, 1.0, MaxwellJuttner(mu=2.0)), NormalizedSpecies(-1.0, 1.0, MaxwellJuttner(mu=2.0)))
    gL = DispersionFunction(pair2, Wavenumber(0.0, 0.85); mode=:L)
    gR = DispersionFunction(pair2, Wavenumber(0.0, 0.85); mode=:R)
    @test gL(0.5 - 0.4im) ≈ gR(0.5 - 0.4im) rtol = 1.0e-10
    # simple sign change across the aperiodic root at ω ≈ −0.479im
    @test real(gL(-0.47im)) * real(gL(-0.485im)) < 0
end
@testitem "Perpendicular factorization: det = O·X, parallel_even gate" begin
    kp = Wavenumber(0.6, 0.0)
    ω = 1.3 + 0.02im
    even = (NormalizedSpecies(1.0, 1.0, Maxwellian(0.05)), NormalizedSpecies(-1.0, 1836.0, Maxwellian(0.05)))
    fO, fX, fD = map(m -> DispersionFunction(even, kp; mode=m), (:O, :X, :det))
    @test fD(ω) ≈ fO(ω) * fX(ω) rtol = 1.0e-12
    # a field-aligned drift recouples the blocks: rejected with opt-in hint
    drift = (NormalizedSpecies(1.0, 1.0, Maxwellian(vth_para=0.05, vd=0.02)), NormalizedSpecies(-1.0, 1836.0, Maxwellian(0.05)))
    @test_throws ArgumentError DispersionFunction(drift, kp; mode=:X)
end

@testitem "Parallel mode solve: ALPS root and close-pass continuation" begin
    using VlasovMaxwellDispersion: successful_retcode
    pair(vdf) = (NormalizedSpecies(1.0, 1.0, vdf), NormalizedSpecies(-1.0, 1.0, vdf))
    p2 = pair(MaxwellJuttner(mu=2.0))
    # ALPS test_relativistic root at exact parallel
    s = solve(DispersionProblem(p2, complex(3.9621e-2, -2.6e-6), Wavenumber(0.0, 0.1); mode=:L))
    @test successful_retcode(s)
    @test real(s.omega) ≈ 0.039185 atol = 1.0e-5
    @test s.resid < 1.0e-10
    # μ=10 aperiodic family through its close pass with the propagating root
    # near k∥ ≈ 3.2–3.3: naive per-k re-seeding jumps branches there; the
    # Continuation predictor gate must not.
    p10 = pair(MaxwellJuttner(mu=10.0))
    ks = [Wavenumber(0.0, kz) for kz in 3.6:-0.1:3.0]
    sc = solve(DispersionProblem(p10, complex(0.0, -2.767), ks; mode=:L))
    @test successful_retcode(sc)
    @test maximum(abs ∘ real, sc.omega) < 1.0e-8         # stays purely imaginary
    @test imag(sc.omega[end]) ≈ -4.334 atol = 5.0e-3     # deep-family value at k=3
end
