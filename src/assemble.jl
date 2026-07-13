# Pointwise perp tensor at node p⊥ before parallel integration.
# M=(q,uq,u²q,p,up); here M=c·Δm. Split into LINEAR FORMS of Δm
# and the Bessel-bilinear assembly (`_In_assemble`), so callers can apply linear
# operations (peel subtraction, Cauchy weights) to 4 scalars instead of 5 moments.
@inline function _In_forms(Δm, px, ω, kz)
    Δ0, Δ1, Δ2, Δ3, Δ4 = Δm
    kzpx = kz * px
    return SA[
        ω * Δ0 - kz * Δ1 + kzpx * Δ3,
        ω * Δ1 - kz * Δ2 + kzpx * Δ4,
        Δ2 - px * Δ4, px * Δ4,
    ]
end
@inline function _In_assemble(F, bvec, nΩ, ω)
    b11, b12, b22, b13, b23, b33 = bvec
    D0, D1 = F[1], F[2]
    zz = b33 * (nΩ * F[3] + ω * F[4])
    return SA[b11 * D0, im * b12 * D0, b13 * D1, b22 * D0, im * b23 * D1, zz]
end
@inline _In_block(Δm, c, bvec, px, ω, kz, nΩ) =
    _In_assemble((2π * c) .* _In_forms(Δm, px, ω, kz), bvec, nΩ, ω)

# Materialize the antisymmetric-paire
@inline _antisymmat(t) =
    @SMatrix [t[1] t[2] t[3]; -t[2] t[4] -t[5]; t[3] t[5] t[6]]

# Symmetric 3×3 from its 6 distinct entries (row-major upper triangle).
@inline _symmat(a11, a12, a13, a22, a23, a33) =
    @SMatrix [a11 a12 a13; a12 a22 a23; a13 a23 a33]

@inline function _perp_Bessel_bilinear(n, a, px)
    z = a * px
    Jm, Jp = besselj(n - 1, z), besselj(n + 1, z)
    b1 = px * (Jm + Jp) / 2
    b2 = px * (Jm - Jp) / 2
    b3 = besselj(n, z)
    return SA[b1 * b1, b1 * b2, b2 * b2, b1 * b3, b2 * b3, b3 * b3]
end

# Fill `out[i]` with the ±nmax ladder of perp bilinear products
function _perp_Bessel_bilinears!(out, a, px)
    z = a * px
    nmax = (length(out) - 1) ÷ 2
    @no_escape begin
        Jv = @alloc(typeof(z), nmax + 2)
        besselj_ladder!(Jv, nmax + 1, z)
        @inbounds for (i, n) in enumerate(-nmax:nmax)
            Jm = _jladder(Jv, n - 1)
            Jp = _jladder(Jv, n + 1)
            Rn = (Jm + Jp) / 2
            Jn = _jladder(Jv, n)
            Jn′ = (Jm - Jp) / 2
            b1 = px * Rn
            b2 = px * Jn′
            out[i] = SA[b1 * b1, b1 * b2, b2 * b2, b1 * Jn, b2 * Jn, Jn * Jn]
        end
    end
    return out
end
