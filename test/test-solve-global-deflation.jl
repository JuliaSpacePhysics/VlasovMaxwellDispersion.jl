@testitem "global solve near ω=0: light-term pole deflated" begin
    # det 𝒟 has a genuine pole at ω=0 (Maxwell light term (kk−k²I)c²/ω², order ≤6).
    # Un-deflated GRPF smears it into spurious "poles" at nonzero ω and drags the
    # neighbouring root estimates off-position; the argument-principle winding of
    # any region containing the origin partially cancels against the nearby roots.
    # The global path must therefore solve det(ω²𝒟) = ω⁶ det𝒟, which is entire.
    mp_me = 1836.15267343
    vA = 1.0e-3
    plasma = (
        NormalizedSpecies(1.0, 1 / vA^2, Maxwellian(vA)),
        NormalizedSpecies(-mp_me, mp_me / vA^2, Maxwellian(vA * sqrt(mp_me))),
    )
    k = Wavenumber(0.0, 0.1)
    region = (-6.0e-4 - 2.0e-4im, 6.0e-4 + 2.0e-4im)

    sol = solve(GlobalDispersionProblem(plasma, k, region), GRPF(; tol = 2.0e-5))

    # Muller-polished references (double Alfvén roots ±kvA and two damped pairs).
    expected = [
        -0.00023608871 - 0.0001811976im,
        -0.00014566586 - 6.2714793e-5im,
        -1.0008e-4 + 0.0im,
        1.0008e-4 + 0.0im,
        0.00014566586 - 6.2714793e-5im,
        0.00023608871 - 0.0001811976im,
    ]
    for t in expected
        @test minimum(abs.(sol.omega .- t)) < 2.0e-5
    end
    # no spurious poles: the only pole in the box is the deflated one at ω=0
    @test isempty(sol.poles)
    # the deflation artifact (zero of order 6−p at exactly ω=0) must be filtered out
    @test all(abs(ω) > 2.0e-5 for ω in sol.omega)
end
