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

struct PeeledQuadGK{T, V}
    lims::T
    ζs::V
end

# per pole subtract g(ζₚ) when it is finite and not too large vs the on-contour scale
# else integrate raw (the far/overflow branch)
function (alg::PeeledQuadGK)(g; side = 1, maxratio = 1.0e6, kw...)
    lo, hi = alg.lims
    ζs = alg.ζs
    gscale = maximum(ζ -> NORM(g(clamp(real(ζ), lo, hi))), ζs)
    gζs = g.(ζs)
    peel = @. all(isfinite, gζs) && NORM(gζs) <= maxratio * gscale
    I = similar(gζs)
    f! = function (y, u)
        gu = g(u)
        return @inbounds for i in eachindex(gζs)
            y[i] = ifelse(peel[i], gu - gζs[i], gu) * inv(u - ζs[i])
        end
    end
    QuadGK.quadgk!(f!, I, lo, hi; norm = v -> maximum(NORM, v), kw...)
    @. I = I + gζs * _lpole_term(ζs, lo, hi, side, peel)
    return I
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

"""
    hilbert(g, ζ, L, U; rtol=1e-9, σ=1)

Landau-causal Cauchy integral `∫_L^U g(v)/(v − ζ) dv` for analytic `g`.

Plemelj split at weakly damped/growing modes to remove singularity:

    ∫_L^U g/(v−ζ) = ∫_L^U (g(v)−g(ζ))/(v−ζ) dv  +  g(ζ)·log((U−ζ)/(L−ζ))  [+ σ·2πi·g(ζ)]

Falls back to the direct integrand when the subtraction is ill-conditioned (see [`_subtract_safe`](@ref)).

`σ = sign(k∥)` orients the contour: the causal (Im ω > 0) side is `σ·Im ζ > 0`, and the
residue `σ·2πi·g(ζ)` is the Landau continuation onto the damped side.
"""
function hilbert(g, ζ, L, U; rtol = 1.0e-9, σ = 1)
    gζ = g(ζ)
    near = _subtract_safe(gζ, abs(g(clamp(real(ζ), L, U))))
    gsub = near ? gζ : zero(gζ)
    reg = QuadGK.quadgk(v -> (g(v) - gsub) / (v - ζ), L, U; rtol)[1]
    return reg + gζ .* _lpole_term(ζ, L, U, σ, near)
end

# Subtracting g(ζ) cancels ~log₁₀(|g(ζ)|/gscale) digits against the analytic log term, and
# g(ζ) overflows outright for strongly damped ζ
@inline _subtract_safe(gζ, gscale) = all(isfinite, gζ) && NORM(gζ) * sqrt(eps(one(gscale))) ≤ gscale

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
