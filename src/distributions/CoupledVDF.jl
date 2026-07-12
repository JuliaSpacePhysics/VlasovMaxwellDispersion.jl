"""
    CoupledVDF(f0; para=(lo,hi), perp=(lo,hi), dgrad=nothing, regime=NonRelativistic())

General analytic gyrotropic VDF `f0(p⊥,p∥)`.

And `para`/`perp` are `(lower, upper)` integration ranges.

`dgrad(p⊥,p∥) -> (∂⊥f0, ∂∥f0)` supplies the gradient and default to autodiff.

Prefer [`SeparableVDF`] when `f0(p⊥,p∥)=f⊥(p⊥)f∥(p∥)`.
"""
struct CoupledVDF{F, Dg, T, R <: Regime} <: AbstractVDF
    f0::F
    dgrad::Dg
    para::Tuple{T, T}
    perp::Tuple{T, T}
    regime::R
end

regime(d::CoupledVDF) = d.regime

@inline _pair(x::Tuple) = x
@inline _pair(x) = (zero(x), x)

function CoupledVDF(f0; para, perp, dgrad = nothing, regime = NonRelativistic())
    plo, phi = promote(float(para[1]), float(para[2]))
    qlo, qhi = oftype(phi, _pair(perp)[1]), oftype(phi, _pair(perp)[2])
    dg = isnothing(dgrad) ? ((q, u) -> _grad2(f0, q, u)) : dgrad
    return CoupledVDF(f0, dg, (plo, phi), (qlo, qhi), regime)
end

function contribution(d::CoupledVDF, s, ω, k; closure = HarmonicSum(), kw...)
    return contribution(prepare(d, closure), s, ω, k; closure, kw...)
end

density(d::CoupledVDF; rtol = 1.0e-9) = 2π * QuadGK.quadgk(
    q -> q * QuadGK.quadgk(u -> d.f0(q, u), d.para...; rtol)[1],
    d.perp...; rtol
)[1]

pperp2_mean(d::CoupledVDF, n = density(d); rtol = 1.0e-3) = 2π * QuadGK.quadgk(
    q -> q^3 * QuadGK.quadgk(u -> d.f0(q, u), d.para...; rtol)[1],
    d.perp...; rtol
)[1] / n

prepare(d::CoupledVDF, closure = HarmonicSum(); kw...) =
    PreparedVDF(d, precompute(regime(d), closure, d; kw...))

precompute(::NonRelativistic, ::Newberger, d; kw...) = (; n = density(d))
function precompute(::NonRelativistic, ::HarmonicSum, d; kw...)
    n = density(d)
    return (; n, pperp2_mean = pperp2_mean(d, n))
end
precompute(::Relativistic, ::Any, d; quad = BoxQuad(_GL24, _GL32), kw...) =
    (; n = density(d), bernstein33 = _bernstein_rel(d, quad))

contribution(c::PreparedVDF, s, ω, k; closure = HarmonicSum(), kw...) =
    _coupled_contribution(closure, regime(c), c, s, complex(float(ω)), k; kw...) / c.cache.n

function _coupled_contribution(::HarmonicSum, ::NonRelativistic, c, s, ω, k; alg = PeeledGK(), norm = NORM, rtol = 1.0e-6)
    d = c.vdf
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    a = kperp / Ω
    nmax = nmax_bessel(a^2 * abs(c.cache.pperp2_mean) / 2)
    ns = (-nmax):nmax
    b2s = similar(ns, SVector{6, typeof(a)})
    X = iszero(kz) ? _coupled_X0(d, ω, Ω, a, ns, b2s; rtol, norm) :
        _coupled_X(alg, d, ω, Ω, kz, a, ns, b2s; rtol, norm)
    return (s.Pi2 / ω^2) * _antisymmat(X)
end

# kz≠0: outer perp quadrature over the ladder primitive.
function _coupled_X(alg, d, ω, Ω, kz, a, ns, b2s; rtol, norm)
    invkz = -1 / kz
    ζs = [(ω - n * Ω) / kz for n in ns]
    ctx = (; lims = d.para, ζs, side = sign(kz), nΩs = ns * Ω, ω, kz)
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

const _GL24 = GaussLegendre(24)
const _GL32 = GaussLegendre(32)

# Relativistic path, sliced in (p⊥,p∥) — docs/relativistic.md.
# Resonance D(p∥) = ωγ − k∥p∥ − nΩ₀ with γ=√(1+p⊥²+p∥²) rationalizes,
#   D·D̃ = A(p∥−p₊)(p∥−p₋),  D̃ = ωγ + k∥p∥ + nΩ₀,  A = ω²−k∥²,
# into two explicit simple poles; the squaring ghost (zero of D̃) carries a null
# residue automatically. Poles cross the real p∥ axis ONLY at Im ω = 0.
# Endpoints |p∥|=P sit where f₀≈0, so no endpoint (rim-type) corrections arise.
# f₀ must be evaluable at complex p∥ (poles sit off-axis for complex ω).
# Validated vs Maxwell–Jüttner (Swanson) to ~1e-5 down to Im ω = −0.15 at μ=2.
function _coupled_contribution(::HarmonicSum, ::Relativistic, c, s, ω, k; quad = BoxQuad(_GL24, _GL32), rtol = 1.0e-6)
    d = c.vdf
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    if imag(ω) < 0 && real(ω)^2 > kz^2
        @warn "damped superluminal ω (|Re ω| > |k∥|): the (p⊥,p∥) integral is not the analytic continuation there (apex branch cut, docs/relativistic.md); evaluate at Im ω ≥ 0 and continue externally" maxlog = 1
    end
    a = kperp / Ω
    qhi = d.perp[2]
    nmax = nmax_bessel(a^2 * qhi^2 / 2)
    f = n -> _coupled_harmonic_rel(n, d, ω, Ω, kz, a, quad)
    X_T = 2π * converge(f; nmax, rtol)
    X = _antisymmat(X_T) .+ _ee33(c.cache.bernstein33)
    return (s.Pi2 / ω^2) * X
end

# Relativistic non-resonant e∥e∥ term without prefactor
# Edge-mapped: γ→q² concentrates nodes near γ=1, p∥→θ half-angle over the resonance ellipse
function _bernstein_rel(d, qs = BoxQuad(_GL24, _GL32))
    γmax = sqrt(1 + max(d.para[1]^2, d.para[2]^2) + d.perp[2]^2)
    acc = quad(qs.outer, 0, 1) do q
        γ = 1 + (γmax - 1) * q^2
        wγ = 2 * (γmax - 1) * q
        umax = sqrt(γ^2 - 1)
        inner = quad(qs.inner, -1, 1) do t
            θ = t * (π / 2)
            u, w = umax .* sincos(θ)
            dpe, dpa = d.dgrad(w, u)
            (π / 2) * ComplexF64(w * u * dpa - u^2 * dpe)
        end
        wγ * inner
    end
    return 2π * acc
end

@inline _ee33(x) = @SMatrix ComplexF64[0 0 0; 0 0 0; 0 0 x]

# Covariant momentum numerator 𝒰 = ω∂_γf+k∥∂_uf at (γ,p∥) with w=p⊥, rewritten via
# ∂_γ|_u=(γ/w)∂_⊥, ∂_u|_γ=∂_∥−(u/w)∂_⊥ ⇒ 𝒰 = k∥∂_∥f + (ωγ−k∥u)/w · ∂_⊥f.
@inline function _U_cov(d, u, w, γ, ω, kz)
    dpe, dpa = d.dgrad(w, u)
    return kz * dpa + dpe * (ω * γ - kz * u) / w
end


const AType = SVector{6, ComplexF64}


# Poles farther than this from the real p∥ segment leave the integrand smooth enough
# for plain adaptive quadrature; nearer (or Landau-crossed) poles are peeled. Also keeps
# the γ=0 artifact roots of D·D̃ (at p∥=±i√(1+p⊥²), |Im|≥1) out of the peeled set.
const _PQ_NEAR = 1.5

# Partial fractions of the rationalized resonance (docs/relativistic.md) at fixed
# m⊥² = 1+p⊥²: D_n·D̃_n = A·u² + B·u + C, so
#   1/D_n = D̃_n·[c₁/(u−p₁) + c₂/(u−p₂)],  c₁₂ = ∓1/√(B²−4AC).
# Vieta gives the second root without cancellation; A→0 sends p₁→∞, marked non-finite
# (its term is O(A); callers drop non-finite poles).
@inline function _Dn_poles(ω, kz, nΩ, m2)
    A = ω^2 - kz^2
    B = -2 * kz * nΩ
    C = ω^2 * m2 - nΩ^2
    sq = sqrt(B^2 - 4 * A * C)
    abs2(B + sq) < abs2(B - sq) && (sq = -sq)
    return (_home_side((-B - sq) / (2A), ω, kz, m2), -1 / sq),
        (_home_side(2 * C / (-B - sq), ω, kz, m2), 1 / sq)
end

# Exactly-real ω leaves a real pole ON the path: a signed zero nudges it to its home
# side (the Im ω→0⁺ limit, slope dp/dω = γ²/(k∥γ−ωp)) so the boundary-value log in
# `_peel_residue` lands on the correct sheet.
@inline function _home_side(p, ω, kz, m2)
    isfinite(p) || return complex(Inf)
    iszero(imag(p)) || return p
    γp = sqrt(complex(m2 + p^2))
    return complex(real(p), sign(real(γp^2 / (kz * γp - ω * p))) * 0.0)
end

# Residue r = c·W(p) of a peeled pole and its analytic across-box term
# r·[log((hi−p)/(lo−p)) + σ·2πi if Landau-crossed]; (0, 0) when not peeled.
# Peel when Landau-crossed (Im ω<0 dragged the pole off its σ-home side)
# or within _PQ_NEAR of the segment (Plemelj subtraction for quadrature health).
#  The squaring ghost peels harmlessly (W(p)=0 ⇒ r=0); γ-artifact roots (γ(p)=0 ⇒ W=∞) are left unpeeled.
@inline function _peel_residue(p, c, W, γof, ν, lo, hi, σ)
    zz = (zero(AType), zero(AType))
    isfinite(p) || return zz
    crossed = ν < 0 && σ * imag(p) < 0 && lo < real(p) < hi
    near = abs(imag(p)) < _PQ_NEAR && lo - _PQ_NEAR < real(p) < hi + _PQ_NEAR
    (crossed || near) || return zz
    γp = γof(p)
    iszero(γp) && return zz   # exact γ-artifact (kz=0, n=0 degenerates both roots to γ(p)=0)
    r = c .* W(p, γp)
    all(isfinite, r) || return zz
    return r, r .* (log((hi - p) / (lo - p)) + (crossed ? σ * 2π * im : 0))
end

# One harmonic of the (p⊥,p∥) box integral: outer Gauss–Legendre in p⊥; per slice
#   ∫ σ𝓣_n/D_n du = ∫ Σᵢ cᵢ·W/(u−pᵢ) du,  W = σ·D̃_n·𝓣_n,  σ = 𝒰·p⊥/γ,
# with peeled poles kept as single fractions (c·W(u)−r)/(u−p) — the split form
# σ𝓣/D − r/(u−p) carries 1/(u−p)² rounding noise near the pole.
function _coupled_harmonic_rel(n, d, ω, Ω, kz, a, qs::BoxQuad)
    plo, phi = d.para
    qlo, qhi = d.perp
    nΩ = n * Ω
    ν = imag(ω)
    σ = sign(kz)
    return quad(qs.outer, qlo, qhi) do q
        m2 = 1 + q^2
        z = a * q
        γof(u) = sqrt(complex(m2 + u^2))
        σof(u, γ) = begin
            dpe, dpa = d.dgrad(q, u)
            (kz * dpa + dpe * (ω * γ - kz * u) / q) * (q / γ)
        end
        Wof = (u, γ) -> (σof(u, γ) * (ω * γ + kz * u + nΩ)) .* _T_n_bare(n, z, u, q)
        (p1, c1), (p2, c2) = _Dn_poles(ω, kz, nΩ, m2)
        r1, lg1 = _peel_residue(p1, c1, Wof, γof, ν, plo, phi, σ)
        r2, lg2 = _peel_residue(p2, c2, Wof, γof, ν, plo, phi, σ)
        reg = quad(qs.inner, plo, phi) do u
            γ = γof(u)
            if isfinite(p1) || isfinite(p2)
                Wu = Wof(u, γ)
                acc = zero(Wu)
                isfinite(p1) && (acc = acc .+ (c1 .* Wu .- r1) ./ (u - p1))
                isfinite(p2) && (acc = acc .+ (c2 .* Wu .- r2) ./ (u - p2))
                acc
            else                       # A=B=0 (ω=±k∥, n=0): quadratic degenerate, no poles
                (σof(u, γ) .* _T_n_bare(n, z, u, q)) ./ (ω * γ - kz * u - nΩ)
            end
        end
        reg .+ lg1 .+ lg2
    end
end

include("qin.jl")
