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

# the exact relativistic fixtures it unlocks (López-artifact investigation,
# experiments/lopez-anomalous-zone/report.typ).
@testitem "Exact relativistic fixtures: ladder law, EM asymptote, band edge" begin
    using SpecialFunctions: besselk
    pair(vdf) = (NormalizedSpecies(1.0, 1.0, vdf), NormalizedSpecies(-1.0, 1.0, vdf))
    for μ in (2.0, 10.0)
        p = pair(MaxwellJuttner(mu=μ))
        # k→0 aperiodic ladder: γ_n = −2μΩ/(π(2n−1)); n=1 member at k=0.05
        # carries only the O(k²) finite-k shift
        s = solve(DispersionProblem(p, complex(0.0, -2μ / π), Wavenumber(0.0, 0.05); mode=:L))
        @test abs(s.omega + im * 2μ / π) < 5.0e-3 * μ
        # superluminal EM branch: ω² − k² → 2Π²K₁(μ)/K₂(μ) on the real axis,
        # where F_L is real beyond the band; located by bisection there. The
        # continued sheet below the axis grows exponentially toward marginal
        # in-band ω (quasimode-ladder interference), so near-real superluminal
        # roots are sought at Im ω = 0, not with damped seeds.
        ωp2 = 2 * besselk(1, μ) / besselk(2, μ)
        fEM = DispersionFunction(p, Wavenumber(0.0, 10.0); mode=:L)
        g(w) = real(fEM(complex(w, 0.0)))
        lo, hi = sqrt(100 + 0.5 * ωp2), sqrt(100 + 1.5 * ωp2)
        @test g(lo) * g(hi) < 0
        for _ in 1:40
            m = (lo + hi) / 2
            g(m) * g(lo) > 0 ? (lo = m) : (hi = m)
        end
        @test ((lo + hi) / 2)^2 - 100 ≈ ωp2 atol = 0.02
    end
    # damped superluminal ω at exactly parallel k: the continued Swanson
    # integral (subluminal-germ sheet, certified against the corrected López
    # continuation, experiments/lopez-anomalous-zone report §4.1). Fixtures are
    # the report's script-09 values.
    p2s = pair(MaxwellJuttner(mu=2.0))
    for (kz, ωref) in (
            (2.5, 2.57742 - 0.46394im),   # A/IC branch beyond the light line
            (3.0, 3.02819 - 0.18609im),   # second family traversing the band
        )
        sS = solve(DispersionProblem(p2s, ωref, Wavenumber(0.0, kz); mode=:L))
        @test abs(sS.omega - ωref) < 2.0e-5
        @test sS.resid < 1.0e-10
    end
    # continuity across the light-line seam (subluminal evaluator ↔
    # continuation): a wrong sheet would jump by O(|F_L|) ~ 0.5 here; the
    # bound sits above the quadrature floor, far below that
    fS = DispersionFunction(p2s, Wavenumber(0.0, 2.0); mode=:L)
    @test abs(fS((2.0 + 1.0e-6) - 0.45im) - fS((2.0 - 1.0e-6) - 0.45im)) < 1.0e-3
    # resonance band edge: for real ω, cyclotron resonance γ_L(ω − k∥v∥) = Ω
    # is solvable iff ω² − k∥² ≤ Ω² (min over |v|<c is √(ω²−k²)); Im F_L
    # must vanish outside the band and not inside
    p2 = pair(MaxwellJuttner(mu=2.0))
    fL = DispersionFunction(p2, Wavenumber(0.0, 2.0); mode=:L)
    @test abs(imag(fL(complex(sqrt(4 + 0.5), 0.0)))) > 1.0e-2   # inside band
    @test abs(imag(fL(complex(sqrt(4 + 1.3), 0.0)))) < 1.0e-6   # outside band
    @test abs(imag(fL(complex(sqrt(4 + 2.0), 0.0)))) < 1.0e-6
end
