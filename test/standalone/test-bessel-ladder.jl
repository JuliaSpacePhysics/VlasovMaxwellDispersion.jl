# besselj_ladder must match SpecialFunctions.besselj to machine precision across
# the WHOLE z range and for every Real eltype. Two regimes meet at √eps(T):
#   z ≥ √eps(T): one Miller downward recurrence (rescaled past floatmax(T))
#   z < √eps(T): leading-term series J_k=(z/2)^k/k!, exact to O(z²) rel — this is
#                where the recurrence's 2n/z would overflow / NaN on denormal z.
using VlasovMaxwellDispersion
using Test
const VMD = VlasovMaxwellDispersion
using .VMD: besselj, besselj_ladder, _jladder

# Type-generic check: tolerances scale with eps(T). Cases span the recurrence band
# and the small-z series band (incl. z below √eps(T) for each type) and z=0.
function check_type(T; atol, rtol)
    cases = [
        (5, 0.3), (8, 2.0), (12, 5.0), (20, 10.0), (30, 25.0), (40, 1.0),
        (50, 0.5), (50, 50.0), (10, 1.0e-3), (3, 0.0),
        (10, sqrt(eps(T)) / 10), (20, sqrt(eps(T)) / 100),   # inside the series band for this T
    ]
    for (M, z) in cases
        zz = T(z)
        v = besselj_ladder(M, zz)
        @test eltype(v) === T
        @test all(isfinite, v)
        for k in 0:M
            @test isapprox(v[k + 1], besselj(k, zz); atol, rtol)
        end
    end
    return
end

@testset "Float64" check_type(Float64; atol = 1.0e-13, rtol = 1.0e-12)
@testset "Float32" check_type(Float32; atol = 1.0f-5, rtol = 1.0f-4)
@testset "BigFloat" check_type(BigFloat; atol = 1.0e-30, rtol = 1.0e-28)

# Denormal / zero / negative z: series branch must stay finite where the bare
# recurrence (2n/z → Inf) would NaN. J_0(0)=1, all higher orders 0.
for (M, z) in [(3, 1.0e-150), (10, 1.0e-300), (3, -1.0e-9), (40, 5.0e-9), (5, 0.0)]
    v = besselj_ladder(M, z)
    @test all(isfinite, v)
    @test v[1] ≈ besselj(0, z) atol = 1.0e-13
    @test maximum(abs(v[k + 1] - besselj(k, z)) for k in 0:M) < 1.0e-13
end

# signed access J_{−k} = (−1)^k J_k
let v = besselj_ladder(8, 3.7)
    for k in -8:8
        @test _jladder(v, k) ≈ besselj(k, 3.7) atol = 1.0e-14
    end
end
