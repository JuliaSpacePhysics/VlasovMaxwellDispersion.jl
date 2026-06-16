# Validate Lerche-Newberger sum-rule: the truncated reference Σ_{|n|≤N} converges
# as N→∞, and slower as z=k̃⊥ρ grows (J_n stays O(1) until |n|≈z).

# `S0 = Σ_n Jₙ²/(a−n) = π J_{−a}J_a / sin(πa)`
# `S1 = Σ_n n Jₙ²/(a−n) = a·S0 − 1`
# `S2 = Σ_n n² Jₙ²/(a−n) = a·S1`
# `SD = Σ_n z Jₙ Jₙ'/(a−n) = (z/2)·π(J_{−a}'J_a + J_{−a}J_a')/sin(πa)`
# `SJp = π J_{−a}'J_a' / sin(πa)` (the `yy`/(Jₙ')² sum equals `SJp + a/z²`)

@testitem "Lerche-Newberger sum rule" begin
    using SpecialFunctions: besselj
    using VlasovMaxwellDispersion: qin_sums

    # direct sums over integer harmonics |n|≤N with Real-order besselj
    function truncated_sums(a, z, N)
        S0 = zero(complex(a));
        S1 = zero(complex(a));
        S2 = zero(complex(a));
        SD = zero(complex(a))
        for n in (-N):N
            Jn = besselj(n, z)
            Jnp = (besselj(n - 1, z) - besselj(n + 1, z)) / 2  # J_n'(z)
            d = a - n
            w = Jn^2 / d
            S0 += w
            S1 += n * w
            S2 += n^2 * w
            SD += z * Jn * Jnp / d
        end
        return S0, S1, S2, SD
    end

    # nonzero Im(a) is the Landau offset that keeps the sum finite (no integer
    # resonance) and exercises the complex-order Bessel path.
    a = 1.37 + 0.21im
    Nbig = 400  # "to-large-N" reference: J_n ≈ 0 for |n| ≫ z here

    @testset "z = $z (k̃⊥ρ)" for z in (1.0, 2.0, 3.0)
        S = qin_sums(a, z)
        Sref = truncated_sums(a, z, Nbig)
        for (c, r) in zip(S, Sref)
            @test isapprox(c, r; rtol=1e-6, atol=1e-10)
        end
    end

    # Truncation error at FIXED N grows with z
    Nfix = 6
    errs = map((1.0, 2.0, 3.0)) do z
        Strunc, = truncated_sums(a, z, Nfix)
        Sexact, = qin_sums(a, z)
        abs(Strunc - Sexact)
    end
    @test errs[1] < errs[2] < errs[3]
end
