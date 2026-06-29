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

# Fused harmonic loop for two analytic factors. The per-harmonic `para_moments`
# (hilbert→quadgk) and `perp_moments` (quadgk) integrands are harmonic-independent
# in `f/f′` — only the resonance pole ζ_n and the Bessel factor change with n. So
# collapse the whole ±nmax ladder into ONE u-quadrature + ONE v-quadrature, each
# vector-valued over n, dropping f/f′ evals from O(nodes·nmax) to O(nodes). Mirrors
# `_coupled_perp`. Fixed `-nmax:nmax` (cap from `nmax_bessel`) replaces `converge`.
function _separable_harmonics(para::AnalyticFactor, perp::AnalyticFactor, β, ω, Ω, kz, rtol)
    ns = (-nmax_harm(perp, β)):nmax_harm(perp, β)
    Ms = _para_moments_all(para, ω, kz, Ω, ns; rtol)
    P∂s, PFs = _perp_moments_all(perp, ns, β; rtol)
    acc = zero(SMatrix{3, 3, ComplexF64})
    @inbounds for i in eachindex(ns)
        acc += _chi_mblock(Ms[i], P∂s[i], PFs[i], ω, kz, ns[i] * Ω)
    end
    return acc
end

@inline _gpar(p, u) = (fp = p.fdf(u); SVector(fp[1], u * fp[1], u^2 * fp[1], fp[2], u * fp[2]))

# All parallel moments M_n in one u-quadrature: the Plemelj split per pole reuses
# the single g(u)=[f,uf,u²f,f′,uf′] evaluation across every harmonic.
function _para_moments_all(p::AnalyticFactor, ω, kz, Ω, ns; rtol = 1.0e-8)
    ζs = [(ω - n * Ω) / kz for n in ns]
    gζs = [_gpar(p, ζ) for ζ in ζs]
    reg = first(
        QuadGK.quadgk(p.lo, p.hi; rtol, norm = x -> maximum(_relsize, x)) do u
            g = _gpar(p, u)
            [(g - gζs[i]) / (u - ζs[i]) for i in eachindex(ns)]
        end
    )
    return [(-1 / kz) .* (reg[i] + gζs[i] .* _landau_logfac(ζs[i], p.lo, p.hi)) for i in eachindex(ns)]
end

# 12 distinct bvec⊗bvec entries (∂F slice, then F slice) at one perp node, given the
# Bessel triplet (J_{n−1},J_n,J_{n+1}); shared f⊥/f⊥′ across n.
@inline function _perp12(v, Jm, Jn, Jp, fq, dfq, vfq)
    b1, b2, b3 = v * (Jm + Jp) / 2, v * (Jm - Jp) / 2, Jn
    k11, k12, k13, k22, k23, k33 = b1^2, b1 * b2, b1 * b3, b2^2, b2 * b3, b3^2
    return SVector(
        dfq * k11, dfq * k12, dfq * k13, dfq * k22, dfq * k23, dfq * k33,
        vfq * k11, vfq * k12, vfq * k13, vfq * k22, vfq * k23, vfq * k33
    )
end

# All perp moments (P∂_n, PF_n) in one v-quadrature; one Bessel ladder per node
# feeds every harmonic instead of three besselj per (n,node).
function _perp_moments_all(p::AnalyticFactor, ns, β; rtol = 1.0e-8)
    kmin = first(ns) - 1
    P = first(
        QuadGK.quadgk(p.lo, p.hi; rtol, norm = x -> maximum(_relsize, x)) do v
            z = β * v
            Jv = [besselj(k, z) for k in kmin:(last(ns) + 1)]   # J_k ladder, shared over n
            fq, dfq = p.fdf(v)
            vfq = v * fq
            [_perp12(v, Jv[n - 1 - kmin + 1], Jv[n - kmin + 1], Jv[n + 1 - kmin + 1], fq, dfq, vfq) for n in ns]
        end
    )
    P∂s = [2π .* _symmat(Pi[1], Pi[2], Pi[3], Pi[4], Pi[5], Pi[6]) for Pi in P]
    PFs = [2π .* _symmat(Pi[7], Pi[8], Pi[9], Pi[10], Pi[11], Pi[12]) for Pi in P]
    return P∂s, PFs
end