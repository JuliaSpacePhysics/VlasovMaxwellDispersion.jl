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

# Fused harmonic loop for two analytic factors. f/f′ are harmonic-independent — only
# the resonance pole ζ_n and the Bessel factor change with n.
function _separable_harmonics(para::AnalyticFactor, perp::AnalyticFactor, β, ω, Ω, kz; rtol, norm = NORM)
    nmax = nmax_harm(perp, β)
    ns = -nmax:nmax
    Ms = _para_moments_all(para, ω, kz, Ω, ns; rtol)
    M = last(ns) + 1
    return @no_escape begin
        Jv = @alloc(typeof(β), M + 1)
        QuadGK.quadgk(perp.lo, perp.hi; rtol, norm) do v
            z = β * v
            fq, dfq = perp.fdf(v)
            vfq = v * fq
            acc = zero(AType)
            besselj_ladder!(Jv, M, z)        # J_0..J_{nmax+1} in one recurrence, signed-indexed
            @inbounds for (i, n) in enumerate(ns)
                Jm, Jn, Jp = _jladder(Jv, n - 1), _jladder(Jv, n), _jladder(Jv, n + 1)
                # bvec=(v⊥Rn, v⊥Jn′, Jn), Rn=½(J_{n−1}+J_{n+1}); K=bvec⊗bvec shared by ∂F/F slices
                b1, b2, b3 = v * (Jm + Jp) / 2, v * (Jm - Jp) / 2, Jn
                K = _symmat(b1^2, b1 * b2, b1 * b3, b2^2, b2 * b3, b3^2)
                acc += _chi_mblock(Ms[i], (2π * dfq) .* K, (2π * vfq) .* K, ω, kz, n * Ω)
            end
            acc
        end[1]
    end
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
    Is = PeeledQuadGK((p.lo, p.hi), ζs)(u -> _gpar(p, u); side = Int(sign(kz)), rtol)
    return (-1 / kz) .* Is
end
