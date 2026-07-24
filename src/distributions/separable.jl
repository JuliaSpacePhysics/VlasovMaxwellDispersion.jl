# Duck-typed contract (no factor supertypes).
# A *parallel* factor defines
#   para_moments(para, Δ, kz)       -> M = (MF0,MF1,MF2,MT0,MT1)
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

parallel_even(d::Separable) = parallel_even(d.fpara)
(d::Separable)(q, u) = d.fperp(q) * d.fpara(u)

# Analytic 1-D `f` over `[lo,hi]`, with `fdf(x)->(f,f′)`
struct AnalyticFactor{F, FD, T}
    f::F
    fdf::FD
    lo::T
    hi::T
end

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

function _plan_perp_moments(perp, β, nmax, rtol)
    return [perp_moments(perp, n, β) for n in (-nmax):nmax]
end

_plan_perp_moments(perp::AnalyticFactor, β, nmax, rtol) =
    _quad_perp_moments(perp.fdf, perp.lo, perp.hi, β, nmax, rtol)

function _quad_perp_moments(fdf, lo, hi, β, nmax, rtol)
    return map((-nmax):nmax) do n
        norm(m) = max(maximum(abs, m[1]), maximum(abs, m[2]))
        m = QuadGK.quadgk(lo, hi; rtol, norm) do v
            f, df = fdf(v)
            b = _perp_Bessel_bilinear(n, β, v)
            K = _symmat(b[1], b[2], b[4], b[3], b[5], b[6])
            SVector((2π * df) .* K, (2π * v * f) .* K)
        end[1]
        return m[1], m[2]
    end
end

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

struct SeparablePlan{C, S, M, Z, R}
    prepared::C
    species::S
    moments::M
    kz::Z
    rtol::R
end

plan_contribution(d::Separable, s, k; kw...) =
    plan_contribution(prepare(d), s, k; kw...)

function plan_contribution(
        c::PreparedVDF{<:Separable}, s, k; rtol = 1.0e-8, kwargs...
    )
    d = c.vdf
    β = perp(k) / s.Omega
    fperp = perp_setup(d.fperp, β)
    nmax = _nmax(c.cache, fperp, β)
    moments = _plan_perp_moments(fperp, β, nmax, rtol)
    return SeparablePlan(c, s, moments, para(k), rtol)
end

function (p::SeparablePlan)(ω)
    c, s, kz = p.prepared, p.species, p.kz
    Ω = s.Omega
    nmax = length(p.moments) ÷ 2
    ns = (-nmax):nmax
    Ms = _para_moments_all(c.vdf.fpara, ω, kz, Ω, ns; rtol = p.rtol)
    X = sum(eachindex(p.moments)) do i
        n = ns[i]
        P∂, PF = p.moments[i]
        _chi_mblock(Ms[i], P∂, PF, ω, kz, n * Ω)
    end
    return _antisymmat((s.Pi2 / (c.cache.n * ω^2)) * X)
end

# Magnetized susceptibility of a separable f = f⊥⊗f∥, summed over cyclotron harmonics.
contribution(c::PreparedVDF{<:Separable}, s, ω, k; kw...) =
    plan_contribution(c, s, k; kw...)(ω)

# Ring/table/Kappa factors supply a closed-form `nmax_harm`; AnalyticFactor uses cached ⟨v⊥²⟩.
_nmax(_cache, fperp, β) = nmax_harm(fperp, β)
_nmax(cache, ::AnalyticFactor, β) = nmax_bessel(β^2 * abs(cache.perp²) / 2)

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
