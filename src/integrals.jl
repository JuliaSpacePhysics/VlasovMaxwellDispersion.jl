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
    hilbert(g, ζ, L, U; rtol=1e-9) -> Complex

Landau-causal Cauchy integral `∫_L^U g(v)/(v − ζ) dv` for analytic `g`.

Plemelj split at weakly damped/growing modes to remove singularity:

    ∫_L^U g/(v−ζ) = ∫_L^U (g(v)−g(ζ))/(v−ζ) dv  +  g(ζ)·log((U−ζ)/(L−ζ))  [+ 2πi·g(ζ)]

Falls back to the direct integrand when the subtraction is ill-conditioned (see [`_subtract_safe`](@ref)).

The residue `2πi·g(ζ)` is the Landau continuation onto the damped side.
"""
function hilbert(g, ζ, L, U; rtol = 1.0e-9)
    gζ = g(ζ)
    near = _subtract_safe(gζ, abs(g(clamp(real(ζ), L, U))))
    gsub = near ? gζ : zero(gζ)
    reg = QuadGK.quadgk(v -> (g(v) - gsub) / (v - ζ), L, U; rtol)[1]
    return reg + _pole_corr(near, gζ, ζ, L, U)
end

# Subtracting g(ζ) cancels ~log₁₀(|g(ζ)|/gscale) digits against the analytic log term, and
# g(ζ) overflows outright for strongly damped ζ
@inline _subtract_safe(gζ, gscale) = all(isfinite, gζ) && _relsize(gζ) * sqrt(eps(one(gscale))) ≤ gscale

@inline function _pole_corr(near, gζ, ζ, lo, hi)
    near && return gζ .* _landau_logfac(ζ, lo, hi)
    return _landau_active(ζ, lo, hi) ? gζ .* (2π * im) : zero(gζ)
end

# Assumes the `kz>0` convention, so `Im ζ<0 ⟺ Im ω<0`
@inline _landau_active(ζ, lo, hi) = imag(ζ) < 0 && lo < real(ζ) < hi

@inline function _landau_logfac(ζ, lo, hi)
    logfac = log((hi - ζ) / (lo - ζ))
    return _landau_active(ζ, lo, hi) ? logfac + 2π * im : logfac
end

function converge(f, nmin::Integer; rtol, nmax::Integer = 200)
    total = f(0)
    n = 1
    while n <= nmax
        shell = f(n) + f(-n)
        total += shell
        if n >= nmin && _relsize(shell) <= rtol * _relsize(total)
            break
        end
        n += 1
    end
    return total
end

converge(f; kw...) = converge(f, 1; kw...)


@inline _relsize(x::Number) = abs(x)
@inline _relsize(x::AbstractArray) = maximum(abs, x)

NORM(x) = maximum(abs, x)

"""
    nmax_bessel(lambda; pad=5) -> Int

Hard harmonic cap from Bessel asymptotics: `J_n` (and `I_n e^{-λ}`) negligible
for `n > b + pad·b^{1/3}` with `b = √(2λ)`.
"""
@inline function nmax_bessel(lambda; pad = 5)
    b = sqrt(2 * lambda)
    return ceil(Int, b + pad * cbrt(b)) + 1
end
