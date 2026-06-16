"""
    CoupledVDF(f0; parlower, parupper, perpupper, dpar=nothing, dperp=nothing)

**Most general** gyrotropic VDF: an arbitrary analytic `f0(p∥,p⊥)`.

`f0` must be evaluable at complex argument (continued onto the Landau contour).

Relativistic species (`regime=Relativistic()`) integrate in `(γ,p∥)`.

Prefer [`SeparableVDF`] when `f0(p∥,p⊥)=f∥(p∥)f⊥(p⊥)`.
"""
struct CoupledVDF{F,Dp,Dq,T} <: AbstractVDF
    f0::F
    dpar::Dp        # ∂f₀/∂p∥
    dperp::Dq       # ∂f₀/∂p⊥
    parlo::T
    parhi::T
    perphi::T
end
function CoupledVDF(
    f0; parlower, parupper, perpupper, dpar=nothing, dperp=nothing, normalize=true
)
    plo, phi = promote(float(parlower), float(parupper))
    qhi = oftype(phi, perpupper)
    n = normalize ?
        2π * QuadGK.quadgk(
        v -> v * QuadGK.quadgk(u -> f0(u, v), plo, phi; rtol=1.0e-9)[1],
        zero(qhi), qhi; rtol=1.0e-9
    )[1] : one(plo)
    fn = (u, v) -> f0(u, v) / n
    dp = isnothing(dpar) ? ((u, v) -> _dwrt1(fn, u, v)) : ((u, v) -> dpar(u, v) / n)
    dq = isnothing(dperp) ? ((u, v) -> _dwrt2(fn, u, v)) : ((u, v) -> dperp(u, v) / n)
    return CoupledVDF(fn, dp, dq, plo, phi, qhi)
end

# Regime trait picks the coordinate:
#   NonRelativistic — (p∥,p⊥), pole ζ=(ω−nΩ)/k∥ fixed; outer ∫dp⊥.
#   Relativistic    — (γ,p∥),  pole p∥=(γω−nΩ)/k∥;     outer ∫dγ.
function contribution(d::CoupledVDF, s::Species, ω, k::Wavenumber; closure::IntegralClosure=HarmonicSum())
    iszero(perp(k)) &&
        throw(ArgumentError("CoupledVDF: magnetized EM tensor needs kperp≠0 (oblique)"))
    return _coupled_contribution(closure, Regime(s), d, s, complex(float(ω)), k)
end

function _coupled_contribution(::HarmonicSum, ::NonRelativistic, d::CoupledVDF, s, ω, k; norm=x -> maximum(abs, x))
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    a = kperp / Ω
    L, U = d.parlo, d.parhi
    p⊥²_mean = 2π * QuadGK.quadgk(
        v -> v^3 * QuadGK.quadgk(u -> d.f0(u, v), L, U; rtol=1.0e-7)[1],
        zero(d.perphi), d.perphi; rtol=1.0e-7
    )[1]
    nmax = nmax_bessel(a^2 * abs(p⊥²_mean) / 2)
    ns = (-nmax):nmax
    χ = first(QuadGK.quadgk(zero(d.perphi), d.perphi; rtol=1.0e-6, norm) do v
        _coupled_perp(v, ns, d, ω, Ω, kz, kperp, a, L, U)
    end)
    return SMatrix{3,3,ComplexF64}((s.Pi2 / ω^2) * χ)
end

# Relativistic (γ,p∥) momentum-space path. Momentum distribution f₀ must be
# evaluable at complex p⊥ (the pole pushes p⊥ off-axis).
# Validated vs Maxwell–Jüttner (Swanson) to ~1e-6 and → bi-Maxwellian as μ→∞.
function _coupled_contribution(::HarmonicSum, ::Relativistic, d::CoupledVDF, s, ω, k)
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    a = kperp / Ω
    γmax = sqrt(1 + max(d.parhi^2, d.parlo^2) + d.perphi^2)
    nmax = nmax_bessel(a^2 * d.perphi^2 / 2)
    f = n -> _coupled_harmonic_rel(n, d, ω, Ω, kz, a, γmax)
    χ = converge(f, 1, 1.0e-6; nmax)
    χ = χ .+ _ee33(_bernstein_rel(d, ω, γmax))
    return SMatrix{3,3,ComplexF64}((s.Pi2 / ω^2) * χ)
end

# Relativistic non-resonant e∥e∥ Bernstein addend 𝒳_B (derivation §5).
# ≡0 for isotropic f₀ (p∥⁻¹∂∥f₀=p⊥⁻¹∂⊥f₀), but O(1) for anisotropic f₀.
# `atol` is load-bearing: for isotropic f₀ the integrand vanishes, so an rtol-only adaptive rule chases relative accuracy on roundoff (~1e-16) and never halts.
function _bernstein_rel(d::CoupledVDF, ω, γmax; rtol=1.0e-7, atol=1.0e-10)
    first(QuadGK.quadgk(one(real(ω)), γmax; rtol, atol) do γ
        umax = sqrt(γ^2 - 1)
        umax > 0 || return zero(ComplexF64)
        2π * first(QuadGK.quadgk(-umax, umax; rtol, atol) do u
            w = sqrt(γ^2 - 1 - u^2)
            ComplexF64(u * d.dpar(u, w) - (u^2 / w) * d.dperp(u, w))
        end)
    end)
end

@inline _ee33(x) = @SMatrix ComplexF64[0 0 0; 0 0 0; 0 0 x]


# Fixed Gauss–Legendre orders for the relativistic path (outer γ, inner p∥). The
# integrand is smooth (single removable pole subtracted), so GL converges fast:
# (48,64) matches Maxwell–Jüttner to ~1e-6 in ~50ms. Very
# sharp/multi-scale f₀ may need higher orders — bump these.
const _GLγ = QuadGK.gauss(48)
const _GLp = QuadGK.gauss(64)

# Covariant momentum numerator 𝒰 = ω∂_γf+k∥∂_uf at (γ,p∥) with w=p⊥, rewritten via
# ∂_γ|_u=(γ/w)∂_⊥, ∂_u|_γ=∂_∥−(u/w)∂_⊥ ⇒ 𝒰 = k∥∂_∥f + (ωγ−k∥u)/w · ∂_⊥f.
@inline _U_cov(d, u, w, γ, ω, kz) = kz * d.dpar(u, w) + d.dperp(u, w) * (ω * γ - kz * u) / w
@inline _rescale(x, w, lo, hi) = ((lo + hi) / 2 + (hi - lo) / 2 * x, (hi - lo) / 2 * w)

# 3×3 relativistic harmonic integrand 2π·𝒰·𝓣_n at (γ,p∥); bare momenta make 𝓣_n
# regular at w=0. Caller passes w=√(γ²−1−u²) (complex off the real p∥ range).
@inline _rel_integrand(u, w, γ, n, a, ω, kz, d) = (2π * _U_cov(d, u, w, γ, ω, kz)) .* _T_n_bare(n, a * w, a, u, w)
@inline _rel_integrand(u, γ, n, a, ω, kz, d) = _rel_integrand(u, sqrt(complex(γ^2 - 1 - u^2)), γ, n, a, ω, kz, d)


# One relativistic harmonic: outer GL over γ∈[1,γmax]; at each γ a single-pole
# parallel Cauchy (Plemelj split + Landau, same invariant as scalar `hilbert`)
# closed by inner GL on the regularized integrand at pole p∥=(γω−nΩ)/k∥.
function _coupled_harmonic_rel(n, d::CoupledVDF, ω, Ω, kz, a, γmax)
    gn, gw = _GLγ
    pn, pw = _GLp
    acc = zero(SMatrix{3,3,ComplexF64})
    for ig in eachindex(gn)
        γ, wγ = _rescale(gn[ig], gw[ig], one(real(ω)), γmax)
        umax = sqrt(γ^2 - 1)
        ζ = (γ * ω - n * Ω) / kz                 # single Landau pole in p∥
        # Peel the pole only when it lands on the real p∥ disk. Off-disk (this n
        # doesn't resonate at this γ) the Plemelj value _rel_integrand(ζ) probes
        # p⊥=√(γ²−1−ζ²) imaginary ⇒ Bessel→Iₙ overflows at large n,k⊥; the
        # integrand is pole-free there, so integrate it directly.
        if -umax < real(ζ) < umax
            nζ = _rel_integrand(ζ, γ, n, a, ω, kz, d)
            reg = zero(nζ)
            for ip in eachindex(pn)
                u, wu = _rescale(pn[ip], pw[ip], -umax, umax)
                reg = reg .+ wu .* ((_rel_integrand(u, γ, n, a, ω, kz, d) .- nζ) ./ (u - ζ))
            end
            inner = reg .+ nζ .* _landau_logfac(ζ, -umax, umax)
        else
            inner = zero(SMatrix{3,3,ComplexF64})
            for ip in eachindex(pn)
                u, wu = _rescale(pn[ip], pw[ip], -umax, umax)
                inner = inner .+ wu .* (_rel_integrand(u, γ, n, a, ω, kz, d) ./ (u - ζ))
            end
        end
        acc = acc .+ wγ .* ((-1 / kz) .* inner)
    end
    return acc
end

# χ(p⊥) for the WHOLE harmonic sum at one perp node
function _coupled_perp(v, ns, d::CoupledVDF, ω, Ω, kz, kperp, a, L, U)
    # Landau–Hilbert for 5 parallel moments: [∂⊥, u·∂⊥, u²·∂⊥, ∂∥, u·∂∥]
    g5(u) = (q=d.dperp(u, v); p=d.dpar(u, v); SVector(q, u * q, u^2 * q, p, u * p))
    ζs = [(ω - n * Ω) / kz for n in ns]
    gζs = [g5(ζ) for ζ in ζs]
    # per-harmonic perp Bessel moments at this p⊥
    ps = map(n -> _perp_bessel_moments(n, a, v), ns)
    # regularized integral part: Σ_n χ_n with the Plemelj removable singularity
    reg = first(QuadGK.quadgk(L, U; rtol=1.0e-7, norm=x -> maximum(abs, x)) do u
        g = g5(u)
        acc = zero(SMatrix{3,3,ComplexF64})
        @inbounds for i in eachindex(ns)
            m = (-1 / kz) .* ((g - gζs[i]) / (u - ζs[i]))
            acc += _chi_mblock((m[1], m[2], m[3], m[4], m[5]), ps[i], ω, kz, kperp, ns[i] / a)
        end
        acc
    end)
    # analytic log-ratio (+ Landau) part, constant in u
    logacc = zero(SMatrix{3,3,ComplexF64})
    @inbounds for i in eachindex(ns)
        m = (-1 / kz) .* (gζs[i] .* _landau_logfac(ζs[i], L, U))
        logacc += _chi_mblock((m[1], m[2], m[3], m[4], m[5]), ps[i], ω, kz, kperp, ns[i] / a)
    end
    return reg + logacc
end

include("qin.jl")