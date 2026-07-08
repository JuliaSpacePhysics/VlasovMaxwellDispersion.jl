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

AnalyticFactor{T}(f, fdf) where {T} = AnalyticFactor(f, fdf, zero(T), T(Inf))

@inline (p::AnalyticFactor)(v) = p.f(v)

function _separable_harmonics(para, perp::AnalyticFactor, β, ω, Ω, kz; kw...)
    return _separable_harmonics_sum_first(para, perp, β, ω, Ω, kz; kw...)
end

@inline function _gpar(p, u)
    fp, dp = p.fdf(u)
    ufp = u * fp
    return SVector(fp, ufp, ufp * u, dp, u * dp)
end

# All parallel moments M_n in one u-quadrature.
function _para_moments_all(p::AnalyticFactor, ω, kz, Ω, ns; rtol = 1.0e-8)
    if iszero(kz)
        # pole-free: one moment integral serves every harmonic, weighted 1/Δ_n
        I = QuadGK.quadgk(u -> _gpar(p, u), p.lo, p.hi; rtol, norm = NORM)[1]
        return [I ./ (ω - n * Ω) for n in ns]
    end
    ζs = [(ω - n * Ω) / kz for n in ns]
    Is = plan_landau((p.lo, p.hi), ζs, sign(kz))(u -> _gpar(p, u); rtol)
    return (-1 / kz) .* Is
end
