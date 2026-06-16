# --- Complex-order Bessel J ---
# SpecialFunctions only ships real order. 
# Hybrid: ascending power series for small |z| (cancellation-safe there) and the divergent large-argument asymptotic
# (truncated at its smallest term) for |z| ≳ 14, where the series loses
# precision to catastrophic cancellation.
const _BESSELJ_ASYM_Z = 14.0

# Ascending series J_ν(z) = Σ_m (−1)^m (z/2)^{2m+ν} / (m! Γ(m+ν+1)).
function _besselj_series(ν, z; maxiters=4000, tol=1.0e-16)
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

# Large-|z| asymptotic (Abramowitz & Stegun 9.2.5). Divergent series ⇒ stop at
# the smallest term. Reliable for |z| ≳ |ν| and |z| ≳ 14.
function _besselj_asym(ν, z; nterms=80)
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

"""
    besselj_complex(ν, z)

Bessel `J_ν(z)` for complex order `ν`. See `_BESSELJ_ASYM_Z` for
the series/asymptotic crossover and its precision limits at large `|z|`.
"""
@inline function besselj_complex(ν::Real, z::Real)
    return isinteger(ν) ? besselj(Int(ν), float(z)) : besselj(float(ν), float(z))
end
@inline function besselj_complex(ν, z)
    return abs(z) >= _BESSELJ_ASYM_Z ? _besselj_asym(complex(ν), complex(z)) :
           _besselj_series(complex(ν), complex(z))
end

"""
    besselj_prime(ν, z) -> J_ν'(z)

Derivative via the standard recurrence `J_ν'(z) = (J_{ν−1}(z) − J_{ν+1}(z))/2`.
Valid for non-integer ν (the case Newberger needs).
"""
@inline besselj_prime(ν, z) = (besselj_complex(ν - 1, z) - besselj_complex(ν + 1, z)) / 2
