# --- Complex-order Bessel J ---
# SpecialFunctions only ships real order.
# Hybrid: ascending power series for small |z| (cancellation-safe there) and the divergent large-argument asymptotic
# (truncated at its smallest term) for |z| ≳ 14, where the series loses
# precision to catastrophic cancellation.
const _BESSELJ_ASYM_Z = 14.0

# Ascending series J_ν(z) = Σ_m (−1)^m (z/2)^{2m+ν} / (m! Γ(m+ν+1)).
function _besselj_series(ν, z; maxiters = 4000, tol = 1.0e-16)
    half = z / 2
    term = half^ν / gamma(ν + 1)
    s = term
    z2 = -half^2
    for m in 1:maxiters
        term *= z2 / (m * (m + ν))
        s += term
        abs(term) <= tol * abs(s) && break
    end
    return s
end

function _besselj_series_deriv(ν, z; maxiters = 4000, tol = 1.0e-16)
    half = z / 2
    term = half^ν / gamma(ν + 1)
    s = term
    sd = ν * term
    z2 = -half^2
    for m in 1:maxiters
        term *= z2 / (m * (m + ν))
        s += term
        sd += (ν + 2m) * term
        abs(term) <= tol * abs(s) && break
    end
    return s, sd / z
end

# Large-|z| asymptotic (Abramowitz & Stegun 9.2.5). Divergent series ⇒ stop at
# the smallest term. Reliable for |z| ≳ |ν| and |z| ≳ 14.
function _besselj_asym(ν, z; nterms = 80)
    μ = 4 * ν^2
    χ = z - (ν / 2 + oftype(real(ν), 1) / 4) * π
    P = one(complex(z))
    Q = zero(complex(z))
    ak = one(complex(z))           # a_0
    minterm = Inf
    for k in 1:nterms
        ak *= (μ - (2k - 1)^2) / (k * 8)
        t = ak / z^k
        at = abs(t)
        at > minterm && break       # past the optimal truncation
        minterm = at
        if iseven(k)
            P += iseven(k ÷ 2) ? t : -t
        else
            Q += iseven((k - 1) ÷ 2) ? t : -t
        end
    end
    return sqrt(2 / (π * z)) * (P * cos(χ) - Q * sin(χ))
end

# Bessel function J_ν(z) and its derivative for complex order ν.
# See `_BESSELJ_ASYM_Z` for the series/asymptotic crossover.
@inline function besselj_deriv(ν, z)
    νc, zc = complex(ν), complex(z)
    return abs(z) < _BESSELJ_ASYM_Z ? _besselj_series_deriv(νc, zc) : _val_dwrt(x -> _besselj_asym(νc, x), zc)
end

function besselj_ladder!(out, M::Integer, z::T) where {T} # out[k+1] = J_k(z), k=0..M
    # Below √eps(T) the one-term series J_k=(z/2)^k/k! is exact to O(z²)≤eps(T) rel,
    # AND it is the safe path on denormal z, where the recurrence's 2n/z → Inf.
    if abs(z) < sqrt(eps(T))
        half = z / 2
        acc = one(T)
        @inbounds out[1] = acc
        @inbounds for k in 1:M
            acc *= half / k
            out[k + 1] = acc
        end
        return out
    end
    # Miller downward recurrence: recurse `J_{n−1}=(2n/z)J_n−J_{n+1}` (stable downward for `n>z`)
    # normalize by the Neumann identity `J_0 + 2(J_2+J_4+…) = 1`.
    # rescale by inv(BIG) whenever it crosses BIG to stay clear of floatmax(T).
    BIG = sqrt(floatmax(T))
    s = inv(BIG)
    base = max(M, ceil(Int, abs(z)))             # seed must clear BOTH the wanted M and the turning point n≈z
    # extra steps decay the seed error in the n>z stable region; digits gained scale
    # with step count, so widen the margin for higher-precision T (≥ the Float64 margin).
    margin = (sqrt(40 * (base + 1)) + 15) * max(1, precision(T) / 53)
    N = base + ceil(Int, margin)
    fkp1 = zero(T)
    fk = sqrt(floatmin(T))                       # tiny seed leaves full exponent range to grow downward
    nrm = zero(T)
    @inbounds for n in N:-1:1
        fkm1 = (2 * n / z) * fk - fkp1
        (n - 1 <= M) && (out[n] = fkm1)
        iseven(n - 1) && (nrm += fkm1)
        fkp1, fk = fk, fkm1
        if abs(fk) > BIG                         # preserve ratios (hence J)
            fk *= s; fkp1 *= s; nrm *= s
            for j in n:(M + 1)
                out[j] *= s
            end
        end
    end
    invn = inv(2 * nrm - out[1])                 # out[1]=J_0 is double-counted in 2·Σ_even
    @inbounds for k in 0:M
        out[k + 1] *= invn
    end
    return out
end


# Whole integer-order ladder `J_0..J_M` at real `z`
# O(M) work vs M+1 independent `besselj` calls
besselj_ladder(M::Integer, z::T) where {T} = besselj_ladder!(Vector{T}(undef, M + 1), M, z)

# Signed access into a `besselj_ladder` result using J_{−k}=(−1)^k J_k.
@inline _jladder(v, k) = k >= 0 ? @inbounds(v[k + 1]) : (iseven(k) ? @inbounds(v[1 - k]) : @inbounds(-v[1 - k]))
