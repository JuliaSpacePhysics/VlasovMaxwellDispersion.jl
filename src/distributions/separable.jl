# Duck-typed contract (no factor supertypes).
# A *parallel* factor defines
#   para_moments(para, ω, kz, nΩ)   -> M = (MF0,MF1,MF2,MT0,MT1)
# A *perpendicular* factor defines
#   perp_setup(perp, β)           -> prepared   (default: itself; rings/tables override)
#   nmax_harm(prepared, β)        -> Int        (harmonic cap from the perp scale)
#   perp_moments(prepared, n, β)  -> (P∂, PF)

"""Separable(fperp, fpar) -> f = fperp ⊗ fpar"""
struct Separable{Q, P} <: AbstractVDF
    fperp::Q
    fpara::P
end

⊗(perp, para) = Separable(perp, para)
(d::Separable)(q, u) = d.fperp(q) * d.fpara(u)

# Analytic 1-D `f` over `[lo,hi]`, with `fdf(x)->(f,f′)`
struct AnalyticFactor{F, FD, T}
    f::F
    fdf::FD
    lo::T
    hi::T
end

AnalyticFactor{T}(f, fdf) where {T} = AnalyticFactor(f, fdf, zero(T), T(Inf))

(p::AnalyticFactor)(v) = p.f(v)

"""
    SeparableVDF(fperp, fpar; para=(lo,hi), perp=(lo,hi), dfpara=nothing, dfperp=nothing)

Arbitrary **separable analytic** VDF `f(p⊥,p∥) = f⊥(p⊥)·f∥(p∥)` for the **full
magnetized EM** susceptibility. Both factors must be
evaluable at complex arguments (continued onto the Landau contour).

`para`/`perp` are `(lower, upper)` integration ranges (a bare `hi` for `perp` means `(0, hi)`).
"""
function SeparableVDF(fperp, fpar; para, perp, dfpara = nothing, dfperp = nothing)
    function _factor(f, df, lo, hi)
        fdf = isnothing(df) ? (x -> _val_dwrt(f, x)) : (x -> (f(x), df(x)))
        return AnalyticFactor(erase_f1(f, hi), erase_fd1(fdf, hi), lo, hi)
    end
    return _factor(fperp, dfperp, _pair(perp)...) ⊗ _factor(fpar, dfpara, para...)
end

# Force generic quadrature path
SeparableVDF(d::Separable; kwargs...) = SeparableVDF(d.fperp, d.fpara; kwargs...)

perp_setup(perp, β) = perp

# Precompute the (ω,k)-independent normalization `n` (∫f=1) and,
# for AnalyticFactor pairs, ⟨v⊥²⟩ setting the Bessel harmonic cap.
# Ring/table factors are self-normalized (n=1) with a closed-form cap from their perp context.
prepare(d::Separable, args...) = PreparedVDF(d, (; n = 1))
function prepare(d::Separable{<:AnalyticFactor, <:AnalyticFactor}, args...; rtol = 1.0e-10)
    p, q = d.fpara, d.fperp
    npara = QuadGK.quadgk(p, p.lo, p.hi; rtol)[1]
    I1, I3 = QuadGK.quadgk(q.lo, q.hi; rtol) do v
        f = q(v)
        SVector(v * f, v^3 * f)
    end[1]
    return PreparedVDF(d, (; n = npara * 2π * I1, perp² = I3 / I1))
end

contribution(d::Separable, s, ω, k; kw...) = contribution(prepare(d), s, ω, k; kw...)

# Magnetized susceptibility of a separable f = f⊥⊗f∥, summed over cyclotron harmonics.
function contribution(c::PreparedVDF{<:Separable}, s, ω, k; rtol = 1.0e-8, kwargs...)
    d = c.vdf
    Ω, kz = s.Omega, para(k)
    β = perp(k) / Ω
    fperp = perp_setup(d.fperp, β)
    X = _separable_harmonics(d.fpara, fperp, β, ω, Ω, kz; rtol, nmax = _nmax(c.cache, fperp, β))
    return _antisymmat((s.Pi2 / (c.cache.n * ω^2)) * X)
end

# Ring/table/Kappa factors supply a closed-form `nmax_harm`; AnalyticFactor uses cached ⟨v⊥²⟩.
_nmax(_cache, fperp, β) = nmax_harm(fperp, β)
_nmax(cache, ::AnalyticFactor, β) = nmax_bessel(β^2 * abs(cache.perp²) / 2)

_separable_harmonics(para, perp::AnalyticFactor, args...; kw...) =
    _separable_harmonics_sum_first(para, perp, args...; kw...)

# Function barrier: `prepared` type is value-dependent
_separable_harmonics(para, perp, args...; kw...) =
    _separable_harmonics_sum_last(para, perp, args...; kw...)

function _separable_harmonics_sum_last(para, perp, β, ω, Ω, kz; rtol, nmax)
    return converge(; nmax, rtol) do n
        nΩ = n * Ω
        Δ = ω - nΩ
        M = para_moments(para, Δ, kz)
        P∂, PF = perp_moments(perp, n, β)
        return _chi_mblock(M, P∂, PF, ω, kz, nΩ)
    end
end

# Fused single-pass harmonic loop: parallel moments Mₙ are v-independent
function _separable_harmonics_sum_first(para, perp, β, ω, Ω, kz; rtol, norm = NORM, nmax)
    ns = -nmax:nmax
    Ms = _para_moments_all(para, ω, kz, Ω, ns; rtol)
    M = last(ns) + 1
    return @no_escape begin
        Jv = @alloc(typeof(β), M + 1)
        QuadGK.quadgk(perp.lo, perp.hi; rtol, norm) do v
            z = β * v
            fq, dfq = perp.fdf(v)
            vfq = v * fq
            besselj_ladder!(Jv, M, z)        # J_0..J_{nmax+1} in one recurrence, signed-indexed
            sum(enumerate(ns)) do (i, n)
                Jm, Jn, Jp = _jladder(Jv, n - 1), _jladder(Jv, n), _jladder(Jv, n + 1)
                # bvec=(v⊥Rn, v⊥Jn′, Jn), Rn=½(J_{n−1}+J_{n+1}); K=bvec⊗bvec shared by ∂F/F slices
                b1, b2, b3 = v * (Jm + Jp) / 2, v * (Jm - Jp) / 2, Jn
                K = _symmat(b1^2, b1 * b2, b1 * b3, b2^2, b2 * b3, b3^2)
                _chi_mblock(Ms[i], (2π * dfq) .* K, (2π * vfq) .* K, ω, kz, n * Ω)
            end
        end[1]
    end
end


_para_moments_all(p, ω, kz, Ω, ns; rtol = 1.0e-8) = map(n -> para_moments(p, ω - n * Ω, kz), ns)

# All parallel moments M_n in one u-quadrature.
function _para_moments_all(p::AnalyticFactor, ω, kz, Ω, ns; rtol = 1.0e-8)
    @inline function _gpar(p, u)
        fp, dp = p.fdf(u)
        ufp = u * fp
        return SVector(fp, ufp, ufp * u, dp, u * dp)
    end

    if iszero(kz)
        # pole-free: one moment integral serves every harmonic, weighted 1/Δ_n
        I = QuadGK.quadgk(u -> _gpar(p, u), p.lo, p.hi; rtol, norm = NORM)[1]
        return [I ./ (ω - n * Ω) for n in ns]
    end
    ζs = [(ω - n * Ω) / kz for n in ns]
    Is = plan_landau((p.lo, p.hi), ζs, sign(kz))(u -> _gpar(p, u); rtol)
    return (-1 / kz) .* Is
end

# Builds one cyclotron-harmonic block χ_n by contracting the perp Bessel tensor
# with the parallel Landau moments (derivation §5.1). Same algebra for every VDF;
# only how the moments are obtained differs (Z/Γ_n closed forms for Maxwellian vs
# `hilbert`+Bessel quadrature for arbitrary f).
#
# The numerator p⊥U splits into a ∂f/∂p⊥ and a ∂f/∂p∥ gradient slice, giving two
# perp Bessel-bilinear matrices and two parallel-moment families:
#   P∂  ← ∫(Bessel)f⊥′    pairs with the f∥ moments M_F^m  (∂⊥ slice)
#   PF  ← ∫(Bessel)f⊥·p⊥  pairs with the f∥′ moments M_T^m  (∂∥ slice)
@inline function _chi_mblock(M, P∂, PF, ω, kz, nΩ)
    MF0, MF1, MF2, MT0, MT1 = M
    # Parallel Landau weights D_m = ω M_F^m − k∥ M_F^{m+1} (∂⊥ slice) and k∥ M_T^m (∂∥ slice).
    # Each tensor entry = (∂⊥ perp bilinear)·wF + (∂∥ perp bilinear)·wT, at order m =
    wF0, wT0 = ω * MF0 - kz * MF1, kz * MT0
    wF1, wT1 = ω * MF1 - kz * MF2, kz * MT1
    xx = P∂[1, 1] * wF0 + PF[1, 1] * wT0
    xy = im * (P∂[1, 2] * wF0 + PF[1, 2] * wT0)
    yy = P∂[2, 2] * wF0 + PF[2, 2] * wT0
    xz = P∂[1, 3] * wF1 + PF[1, 3] * wT1
    zy = im * (P∂[2, 3] * wF1 + PF[2, 3] * wT1)
    zz = nΩ * P∂[3, 3] * MF2 + (ω - nΩ) * PF[3, 3] * MT1   # + non-resonant term
    return SA[xx, xy, xz, yy, zy, zz]
end
