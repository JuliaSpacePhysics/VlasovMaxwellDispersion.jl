NORM(x) = maximum(abs, x)
const _2πim = 2π * im

function _lpole_term(ζ, lo, hi, side, peeled)
    inside = lo < real(ζ) < hi
    cross = (inside && side * imag(ζ) < 0) ? side * _2πim : zero(_2πim)
    peeled || return cross
    return if inside && iszero(imag(ζ))
        complex(log((hi - real(ζ)) / (real(ζ) - lo)), side * π)
    else
        log((hi - ζ) / (lo - ζ)) + cross
    end
end

# ---- 1-D quadrature schemes
abstract type QuadScheme end

struct GaussLegendre{NW} <: QuadScheme
    nw::NW                       # (nodes, weights) === QuadGK.gauss(n)
end
GaussLegendre(n::Integer) = GaussLegendre(QuadGK.gauss(n))

struct AdaptiveGK{K <: NamedTuple} <: QuadScheme
    kw::K
end
AdaptiveGK(; kw...) = AdaptiveGK(NamedTuple(kw))

@inline function quad(f, s::GaussLegendre, lo, hi)
    n, w = s.nw
    mid, half = (lo + hi) / 2, (hi - lo) / 2
    acc = (half * w[1]) * f(mid + half * n[1])
    @inbounds for i in 2:length(n)
        acc = acc + (half * w[i]) * f(mid + half * n[i])
    end
    return acc
end
quad(f, s::AdaptiveGK, lo, hi) = QuadGK.quadgk(f, lo, hi; s.kw...)[1]

# Two composable 1-D schemes for a 2-D box integral (outer × inner).
struct BoxQuad{O <: QuadScheme, I <: QuadScheme}
    outer::O
    inner::I
end


# ---- Landau-causal Cauchy integral
# ∫ g(v)/(v−ζ) dv with the Landau prescription
# σ = sign(k∥) orients the contour: the causal (Im ω > 0) side is `σ·Im ζ > 0`
# residue `σ·2πi·g(ζ)` is the Landau continuation onto the damped side.
abstract type LandauAlg end

"""
Adaptive QuadGK with per-pole subtraction for poles within `_PEEL_BAND·(hi−lo)` of the
real axis, which removes the near-singularity for weakly damped/growing modes:

    ∫_L^U g/(v−ζ) = ∫_L^U (g(v)−g(ζ))/(v−ζ) dv  +  g(ζ)·log((U−ζ)/(L−ζ))  [+ σ·2πi·g(ζ)]

Farther poles and subtraction that is ill-conditioned within the band integrate directly.
"""
struct PeeledGK <: LandauAlg end

# Peel band: |Im ζ| ≤ _PEEL_BAND·(hi−lo). Beyond it the direct integrand is smooth at
# quadrature scale (≲ log₂(1/_PEEL_BAND) extra bisection levels).
# Cf. _PQ_NEAR for relativistic case.
const _PEEL_BAND = 1 / 64

# `alg` selects the numerical method — extensible.
struct LandauPlan{A, T, V, S}
    alg::A
    lims::T
    ζs::V
    side::S
end
plan_landau(lims, ζs, side = 1; alg = PeeledGK()) = LandauPlan(alg, lims, ζs, side)

(p::LandauPlan)(g; kw...) = _landau(p.alg, g, p.lims, p.ζs, p.side; kw...)

# Peel only while ≲2 digits cancel between g(ζ) and the analytic log term
@inline _peel(gζ, gscale) = all(isfinite, gζ) && NORM(gζ) ≤ 1.0e2 * gscale

# Nodes within `_NODE_BAND·U` of ζ are EXCLUDED from the fixed-node Cauchy kernel
# wₗ = uwₗ/(uₗ−ζ). Their terms in Σₗwₗφₗ and in φ(ζ)·Σₗwₗ both reach ~uw/|δ| and must cancel
# back to O(1), so keeping them loses eps·uw/|δ| — 2% of A₀ at |δ|=1e-16. The caller
# substitutes the removable limit uwₗ·φ′(ζ) instead (excluded weight returned as `w0`), 
# The band balances the two: below it cancellation dominates, above it the φ″·|δ| truncation does.
const _NODE_BAND = 1.0e-8

# → (Σₗwₗ over kept nodes, Σuwₗ over excluded nodes, any excluded)
function _node_kernel!(w, un, uw, ζ, near)
    W = zero(ζ)
    w0 = zero(eltype(uw))
    hit = false
    @inbounds for l in eachindex(w)
        δ = un[l] - ζ
        if abs(δ) < near
            w[l] = zero(ζ); w0 += uw[l]; hit = true
        else
            w[l] = uw[l] * safe_inv(δ)
            W += w[l]
        end
    end
    return W, w0, hit
end

# Initial QuadGK segments split at the Landau-pole real parts inside (lo,hi)
function _quadgk_pole_segments(ζs, lo, hi)
    bnds = sort!(unique!(push!(clamp.(real.(ζs), lo, hi), lo, hi)))
    return QuadGK.Segment.(@view(bnds[1:(end - 1)]), @view(bnds[2:end]))
end

# Per-pole reduction. Returns `(gsub, gres)`
@inline function _peel_pole(g, ζ, lo, hi, side)
    if abs(imag(ζ)) > (hi - lo) * _PEEL_BAND
        crossed = lo < real(ζ) < hi && side * imag(ζ) < 0
        gζ = crossed ? g(ζ) : zero(g(clamp(real(ζ), lo, hi))) .* one(ζ)
        return zero(gζ), gζ .* (crossed ? side * _2πim : zero(_2πim))
    end
    gζ = g(ζ)
    peeled = _peel(gζ, NORM(g(clamp(real(ζ), lo, hi))))
    lpole = _lpole_term(ζ, lo, hi, side, peeled)
    gres = iszero(lpole) ? zero(gζ) .* _2πim : gζ .* lpole
    return (peeled ? gζ : zero(gζ)), gres
end

function _landau(::PeeledGK, g, lims, ζ::Number, side; kw...)
    lo, hi = lims
    gsub, gres = _peel_pole(g, ζ, lo, hi, side)
    segs = _quadgk_pole_segments([ζ], lo, hi)
    reg = QuadGK.quadgk(v -> (g(v) - gsub) / (v - ζ), lo, hi; eval_segbuf = segs, kw...)[1]
    return reg + gres
end

# conj/abs2 over Base's overflow-safe inv(::ComplexF64)
safe_inv(x) = conj(x) / abs2(x)

function _landau(::PeeledGK, g, lims, ζs::AbstractVector, side; kw...)
    lo, hi = lims
    gsub_res = _peel_pole.(g, ζs, lo, hi, side)
    gsubs = first.(gsub_res)
    I = similar(ζs, typeof(first(gsub_res)[2]))
    f! = function (y, u)
        gu = g(u)
        return @inbounds for i in eachindex(ζs)
            y[i] = (gu - gsubs[i]) * safe_inv(u - ζs[i])
        end
    end
    segs = _quadgk_pole_segments(ζs, lo, hi)
    QuadGK.quadgk!(f!, I, lo, hi; eval_segbuf = segs, norm = v -> maximum(NORM, v), kw...)
    @inbounds for i in eachindex(I, gsub_res)
        I[i] += last(gsub_res[i])
    end
    return I
end

# ---- Parallel-stage primitive: the harmonic-ladder Landau reduction.
# At one perp node p⊥=v with Bessel bilinears `b2s`, a `LandauAlg` computes
#   X(v) = Σᵢ _In_assemble(ℒ[_In_forms∘g](ζᵢ), b2s[i], nᵢΩ, ω),
# with ℒ[h](ζ) = ∫ h(u)/(u−ζ) du. This is the backend seam — a `plan_ladder`
# method plus a call overload per alg. PeeledGK fuses the pole sum into one
# adaptive peeled integrand. The plan holds the per-(ω,k) context and all
# buffers, so calling it per perp node is allocation-free. Default backend:
# `alg` kwarg upstream.
struct LadderPlan{A, C, W}
    alg::A
    ctx::C      # (; lims, ζs, side, nΩs, ω, kz)
    ws::W       # alg-specific workspace/precomputation
end

function plan_ladder(alg::PeeledGK, ctx; rtol)
    CT = eltype(ctx.ζs)
    Fs = Vector{SVector{4, CT}}(undef, length(ctx.ζs))
    segbuf = QuadGK.alloc_segbuf(typeof(float(ctx.lims[1])), SVector{6, CT}, real(CT))
    edges = _quadgk_pole_segments(ctx.ζs, ctx.lims...)
    return LadderPlan(alg, ctx, (; Fs, segbuf, rtol, edges))
end

function (p::LadderPlan{PeeledGK})(g::G, v, b2s) where {G}
    (; lims, ζs, side, nΩs, ω, kz) = p.ctx
    (; Fs, segbuf, rtol, edges) = p.ws
    lo, hi = lims
    X0 = zero(SVector{6, eltype(eltype(Fs))})
    Xan = X0
    @inbounds for i in eachindex(ζs)
        ζ = ζs[i]
        gsub, gres = _peel_pole(g, ζ, lo, hi, side)
        Fs[i] = _In_forms(gsub, v, ω, kz)   # peel subtrahend, stored as its forms
        lp = _In_forms(gres, v, ω, kz)
        Xan = Xan .+ _In_assemble(lp, b2s[i], nΩs[i], ω)
    end
    reg = QuadGK.quadgk(lo, hi; rtol, norm = NORM, segbuf, eval_segbuf = edges) do u
        Fu = _In_forms(g(u), v, ω, kz)
        acc = X0
        @inbounds for i in eachindex(ζs)
            F = (Fu .- Fs[i]) .* safe_inv(u - ζs[i])
            acc = acc .+ _In_assemble(F, b2s[i], nΩs[i], ω)
        end
        acc
    end[1]
    return reg .+ Xan
end

# Scaled Bessel moment `Γ_n(λ) = I_n(λ) e^{-λ}` from perp gyro-averaging.
# `λ = (k⊥ v_th⊥ / Ω_s)^2 / 2`. Uses scaled modified Bessel `besselix`.
@inline Gamma_n(n, lambda) = besselix(n, lambda)

# Precomputed Γ_k(λ) table with signed indexing `Γ[k]=Γ_{|k|}(λ)` (since I_k=I_{-k}). Built
# once per perp setup so the harmonic loop reuses besselix values across ±n shells instead of
# recomputing them. `kmax` must cover the largest |k| any harmonic reaches.
struct GammaTable{T}
    v::Vector{T}    # v[i] = Gamma_n(i-1, λ)
end
GammaTable(λ, kmax::Integer) = GammaTable([Gamma_n(k, λ) for k in 0:kmax])
@inline Base.getindex(t::GammaTable, k::Integer) = @inbounds t.v[abs(k) + 1]

"""Sugar for `plan_landau((L,U), ζ; σ)(g; rtol)`"""
landau(g, ζ, L, U; side = 1, kw...) = plan_landau((L, U), ζ, side)(g; kw...)

function converge(f, nmin::Integer; rtol, nmax::Integer = 200)
    total = f(0)
    n = 1
    while n <= nmax
        shell = f(n) + f(-n)
        total += shell
        if n >= nmin && NORM(shell) <= rtol * NORM(total)
            break
        end
        n += 1
    end
    return total
end

converge(f; kw...) = converge(f, 1; kw...)

"""
    nmax_bessel(lambda; pad=5) -> Int

Hard harmonic cap from Bessel asymptotics: `J_n` negligible
for `n > b + pad·b^{1/3}` with `b = √(2λ)`.
"""
@inline function nmax_bessel(lambda; pad = 5)
    b = sqrt(2 * lambda)
    return ceil(Int, b + pad * cbrt(b)) + 1
end
