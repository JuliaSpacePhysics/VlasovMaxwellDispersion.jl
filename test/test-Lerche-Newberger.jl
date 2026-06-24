# Validate Lerche-Newberger sum-rule: the truncated reference Σ_{|n|≤N} converges
# as N→∞, and slower as z=k̃⊥ρ grows (J_n stays O(1) until |n|≈z).
#
# `qin_sigmas` returns the regularized quartet the kernel needs, each entire in z:
# `σ0 = Σ_n Jₙ²/(a−n)            = π J_{−a}J_a/sin(πa)`
# `σ1 = (Σ_n n Jₙ²/(a−n))/z²     = (a·σ0 − 1)/z²`
# `σD = (Σ_n z Jₙ Jₙ'/(a−n))/z²  = (z/2)π(J_{−a}'J_a + J_{−a}J_a')/sin(πa) / z²`
# `σJ = Σ_n (Jₙ')²/(a−n)         = π J_{−a}'J_a'/sin(πa) + a/z²`
# Covers both branches: z≥1 closed Bessel, z<1 the entire z²-series.

@testitem "Lerche-Newberger sum rule" begin
    using SpecialFunctions: besselj
    using VlasovMaxwellDispersion: qin_sigmas

    # regularized quartet (σ0,σ1,σD,σJ) from direct sums over integer harmonics |n|≤N
    # using real-order besselj
    function truncated_sums(a, z, N)
        σ0 = σ1 = σD = σJ = zero(complex(a))
        for n in (-N):N
            Jn = besselj(n, z)
            Jnp = (besselj(n - 1, z) - besselj(n + 1, z)) / 2  # J_n'(z)
            d = a - n
            w = Jn^2 / d
            σ0 += w
            σ1 += n * w
            σD += z * Jn * Jnp / d
            σJ += Jnp^2 / d
        end
        z2 = z^2
        return σ0, σ1 / z2, σD / z2, σJ
    end

    # nonzero Im(a) is the Landau offset that keeps the sum finite (no integer
    # resonance) and exercises the complex-order Bessel path.
    a = 1.37 + 0.21im
    Nbig = 400

    @testset "z = $z (k̃⊥ρ)" for z in (0.3, 0.7, 1.0, 2.0, 3.0)
        S = qin_sigmas(a, z)
        Sref = truncated_sums(a, z, Nbig)
        for (c, r) in zip(S, Sref)
            @test isapprox(c, r; rtol = 1.0e-6, atol = 1.0e-10)
        end
    end

    Nfix = 6
    errs = map((1.0, 2.0, 3.0)) do z
        Strunc = truncated_sums(a, z, Nfix)[1]
        Sexact = qin_sigmas(a, z)[1]
        abs(Strunc - Sexact)
    end
    @test errs[1] < errs[2] < errs[3]
end
