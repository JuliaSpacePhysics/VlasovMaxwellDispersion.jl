"""
    SeparableVDF(fperp, fpar; para=(lo,hi), perp=(lo,hi), dfpara=nothing, dfperp=nothing)
    SeparableVDF(d::Separable; para, perp)

Arbitrary **separable analytic** VDF `f(p⊥,p∥) = f⊥(p⊥)·f∥(p∥)` for the **full
magnetized EM** susceptibility. Both factors must be
evaluable at complex arguments (continued onto the Landau contour).

`para`/`perp` are `(lower, upper)` integration ranges (a bare `hi` for `perp` means `(0, hi)`).
Each factor is wrapped as an [`AnalyticFactor`] and fed through [`Separable`] harmonic loop.
"""
function SeparableVDF(
        fperp, fpar; para, perp,
        dfpara = nothing, dfperp = nothing, normalize = true
    )
    plo, phi = promote(float(para[1]), float(para[2]))
    qlo, qhi = oftype(phi, _pair(perp)[1]), oftype(phi, _pair(perp)[2])
    np = normalize ? QuadGK.quadgk(fpar, plo, phi; rtol = 1.0e-10)[1] : one(plo)
    nq = normalize ? 2π * QuadGK.quadgk(v -> fperp(v) * v, qlo, qhi; rtol = 1.0e-10)[1] : one(plo)
    fp = u -> fpar(u) / np
    fq = v -> fperp(v) / nq
    return Separable(_factor(fq, dfperp, nq, qlo, qhi), _factor(fp, dfpara, np, plo, phi))
end

function _factor(f, df, n, lo, hi)
    fdf = isnothing(df) ? (x -> _val_dwrt(f, x)) : (x -> (f(x), df(x) / n))
    return AnalyticFactor(f, fdf, lo, hi)
end

# Force generic quadrature path
SeparableVDF(d::Separable; kwargs...) = SeparableVDF(d.fperp, d.fpara; kwargs...)

# Generic separable factor: an arbitrary analytic 1-D `f` over `[lo,hi]`, with `fdf(x)->(f,f′)`
struct AnalyticFactor{F, FD, T}
    f::F
    fdf::FD
    lo::T
    hi::T
end
@inline (p::AnalyticFactor)(v) = p.f(v)

function para_moments(p::AnalyticFactor, ω, kz, nΩ)
    ζ = (ω - nΩ) / kz
    return (-1 / kz) .* hilbert(ζ; lower = p.lo, upper = p.hi) do u
        fp, dp = p.fdf(u)
        SVector(fp, u * fp, u^2 * fp, dp, u * dp)
    end
end

# 6 distinct entries of bvec⊗bvec in each of the f⊥′ (∂F) and v⊥·f⊥ (F) slices 
# where bvec=(v⊥Rn, v⊥Jn′, Jn) and Rn=½(J_{n−1}+J_{n+1}).
function perp_moments(p, n, β; rtol = 1.0e-8)
    function perptri(v)
        z = β * v
        Jm, Jp1 = besselj(n - 1, z), besselj(n + 1, z)
        bvec = SVector(v * (Jm + Jp1) / 2, v * (Jm - Jp1) / 2, besselj(n, z))
        k11, k12, k13, k22, k23, k33 =
            bvec[1]^2, bvec[1] * bvec[2], bvec[1] * bvec[3], bvec[2]^2, bvec[2] * bvec[3], bvec[3]^2
        fq, dfq = p.fdf(v); vfq = v * fq
        SVector(
            dfq * k11, dfq * k12, dfq * k13, dfq * k22, dfq * k23, dfq * k33,
            vfq * k11, vfq * k12, vfq * k13, vfq * k22, vfq * k23, vfq * k33
        )
    end
    P = 2π .* QuadGK.quadgk(perptri, p.lo, p.hi; rtol)[1]
    return _symmat(P[1], P[2], P[3], P[4], P[5], P[6]), _symmat(P[7], P[8], P[9], P[10], P[11], P[12])
end