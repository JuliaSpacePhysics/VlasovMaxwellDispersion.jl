"""
    CoupledVDF(f0; para, perp, dgrad=nothing, coords=:momentum, regime=NonRelativistic())

Gyrotropic VDF `f0(p⊥, p∥)`, integrated over the box 
`p⊥ ∈ perp`, `p∥ ∈ para` (each a `(lo, hi)` tuple).

`coords` names the variables of `f0` and its gradient `dgrad`:

- `:momentum` (default): `f0(p⊥, p∥)`, `dgrad(p⊥, p∥) = (∂_⊥f, ∂_∥f)`.
- `:energy` (relativistic): `f0(γ, p∥)`, `dgrad(γ, p∥) = (∂_γf, ∂_∥f)`, with the Lorentz factor `γ = √(1 + p⊥² + p∥²)`.

`dgrad` defaults to autodiff of `f0`.

Prefer [`SeparableVDF`] when `f(p⊥, p∥) = f⊥(p⊥) f∥(p∥)`.
"""
struct CoupledVDF{F,Dg,De,T,R<:Regime} <: AbstractVDF
    f0::F
    dgrad::Dg      # (∂_⊥f, ∂_∥f)
    denergy::De    # (∂_γf, ∂_∥f)
    para::Tuple{T,T}
    perp::Tuple{T,T}
    regime::R
end

regime(d::CoupledVDF) = d.regime

@inline _pair(x::Tuple) = x
@inline _pair(x) = (zero(x), x)

#   ∂_⊥f = (p⊥/γ) ∂_γf,   ∂_∥f = ∂_∥f|_γ + (p∥/γ) ∂_γf.
@inline function _mom_from_energy(denergy, q, u)
    γ = sqrt(1 + q^2 + u^2)
    dγ, du = denergy(γ, u)
    return ((q / γ) * dγ, du + (u / γ) * dγ)
end

function CoupledVDF(f0; para, perp, dgrad=nothing, coords=:momentum, regime=NonRelativistic())
    plo, phi = promote(float(para[1]), float(para[2]))
    qlo, qhi = oftype(phi, _pair(perp)[1]), oftype(phi, _pair(perp)[2])
    if coords === :energy
        denergy = @something dgrad (γ, u) -> _grad2(f0, γ, u)
        f0mom = (q, u) -> f0(sqrt(1 + q^2 + u^2), u)
        dgmom = (q, u) -> _mom_from_energy(denergy, q, u)
    elseif coords === :momentum
        denergy = nothing
        f0mom = f0
        dgmom = @something dgrad (q, u) -> _grad2(f0, q, u)
    else
        throw(ArgumentError("coords must be :momentum or :energy, got $coords"))
    end
    return CoupledVDF(erase_f2(f0mom, phi), erase_g2(dgmom, phi), denergy, (plo, phi), (qlo, qhi), regime)
end

function contribution(d::CoupledVDF, s, ω, k; closure=HarmonicSum(), kw...)
    return contribution(prepare(d, closure), s, ω, k; closure, kw...)
end

density(d::CoupledVDF; rtol=1.0e-9) = 2π * QuadGK.quadgk(
    q -> q * QuadGK.quadgk(u -> d.f0(q, u), d.para...; rtol)[1],
    d.perp...; rtol
)[1]

pperp2_mean(d::CoupledVDF, n=density(d); rtol=1.0e-3) = 2π * QuadGK.quadgk(
    q -> q^3 * QuadGK.quadgk(u -> d.f0(q, u), d.para...; rtol)[1],
    d.perp...; rtol
)[1] / n

prepare(d::CoupledVDF, closure=HarmonicSum(); kw...) =
    PreparedVDF(d, precompute(regime(d), closure, d; kw...))

precompute(::NonRelativistic, ::Newberger, d; kw...) = (; n=density(d))
function precompute(::NonRelativistic, ::HarmonicSum, d; kw...)
    n = density(d)
    return (; n, pperp2_mean=pperp2_mean(d, n))
end
precompute(::Relativistic, ::Any, d; quad=BoxQuad(_GL24, _GL32), kw...) =
    (; n=density(d), bernstein33=_bernstein_rel(d, quad))

contribution(c::PreparedVDF{<:CoupledVDF}, s, ω, k; closure=HarmonicSum(), kw...) =
    _coupled_contribution(closure, regime(c), c, s, ω, k; kw...) / c.cache.n

function _coupled_contribution(::HarmonicSum, ::NonRelativistic, c, s, ω, k; alg=PeeledGK(), norm=NORM, rtol=1.0e-6)
    d = c.vdf
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    a = kperp / Ω
    nmax = nmax_bessel(a^2 * abs(c.cache.pperp2_mean) / 2)
    ns = (-nmax):nmax
    b2s = similar(ns, SVector{6,typeof(a)})
    X = iszero(kz) ? _coupled_X0(d, ω, Ω, a, ns, b2s; rtol, norm) :
        _coupled_X(alg, d, ω, Ω, kz, a, ns, b2s; rtol, norm)
    return (s.Pi2 / ω^2) * _antisymmat(X)
end

# kz≠0: outer perp quadrature over the ladder primitive.
function _coupled_X(alg, d, ω, Ω, kz, a, ns, b2s; rtol, norm)
    invkz = -1 / kz
    ζs = [(ω - n * Ω) / kz for n in ns]
    ctx = (; lims=d.para, ζs, side=sign(kz), nΩs=ns * Ω, ω, kz)
    plan = plan_ladder(alg, ctx; rtol)
    return QuadGK.quadgk(d.perp...; rtol, norm) do v
        _perp_Bessel_bilinears!(b2s, a, v)
        (2π * invkz) * plan(v, b2s) do u
            q, p = d.dgrad(v, u)
            SVector(q, u * q, u^2 * q, p, u * p)
        end
    end[1]
end

# kz=0: I is harmonic-independent, weight per n by 1/Δ_n = 1/(ω−nΩ)
function _coupled_X0(d, ω, Ω, a, ns, b2s; rtol, norm)
    return QuadGK.quadgk(d.perp...; rtol, norm) do v
        _perp_Bessel_bilinears!(b2s, a, v)
        I = QuadGK.quadgk(d.para...; norm, rtol) do u
            q, p = d.dgrad(v, u)
            SVector(q, u * q, u^2 * q, p, u * p)
        end[1]
        sum(enumerate(ns)) do (i, n)
            _In_block(I, 1 / (ω - n * Ω), b2s[i], v, ω, zero(a), n * Ω)
        end
    end[1]
end

# Relativistic harmonic sum, sliced in (p⊥, p∥) — docs/src/relativistic.typ.
# Resonance D_n(p∥) = ωγ − k∥p∥ − nΩ, γ = √(1+p⊥²+p∥²), rationalizes to two simple
# poles: D_n·D̃_n = A(p∥−p₊)(p∥−p₋), D̃_n = ωγ + k∥p∥ + nΩ, A = ω²−k∥² (the D̃_n zero is
# a ghost carrying null residue). Poles reach the real p∥ axis only at Im ω = 0, and
# |p∥| = P endpoints sit where f₀ ≈ 0 (no rim terms). f₀ must be evaluable at complex
# p∥. Validated vs Maxwell–Jüttner (Swanson) to ~1e-5 down to Im ω = −0.15 at μ = 2.
#
# `path` picks the damped continuation (docs/src/relativistic.typ):
#   :landau — straight box + classic Landau residues; holomorphic per ω half-plane.
#   :cycles — box (Landau off) + transported residue cycles; needs coords = :energy.
#   :auto   — :cycles for damped-superluminal ω when an energy form is available, else
#             :landau (with a wrong-sheet warning at superluminal ω).
function _coupled_contribution(::HarmonicSum, ::Relativistic, c, s, ω, k;
    path=:auto, quad=BoxQuad(_GL24, _GL32), rtol=1.0e-6, scaledUcov=nothing)
    d = c.vdf
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    _use_cycles(path, d, ω, kz) && return _cycle_contribution(c, s, ω, k; quad, rtol, scaledUcov)
    path === :auto && _warn_damped_superluminal(ω, kz)
    a = kperp / Ω
    nmax = nmax_bessel(a^2 * d.perp[2]^2 / 2)
    X_T = 2π * converge(n -> _harmonic_rel(n, d, ω, Ω, kz, a, quad); nmax, rtol)
    X = _antisymmat(X_T) .+ _ee33(c.cache.bernstein33)
    return (s.Pi2 / ω^2) * X
end

# damped-superluminal ⇔ Im ω < 0 and |Re ω| > |k∥| (k∥ ≠ 0)
@inline function _use_cycles(path, d, ω, kz)
    path === :landau && return false
    damped_super = imag(ω) < 0 && real(ω)^2 > kz^2 && !iszero(kz)
    (path === :cycles || (path === :auto && damped_super)) || return false
    if isnothing(d.denergy)
        path === :cycles && throw(ArgumentError("path = :cycles needs coords = :energy (an analytic denergy)"))
        return false
    end
    return true
end

include("qin.jl")
