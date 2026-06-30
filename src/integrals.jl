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
    hilbert(g, ζ; lower, upper, rtol=1e-9) -> Complex

Landau-causal Cauchy integral `∫_lower^upper g(v)/(v − ζ) dv` for any callable
`g` that is *analytic* (evaluable at complex argument). Plemelj split, exact for
all complex `ζ`:

    ∫ g/(v−ζ) = ∫ (g(v)−g(ζ))/(v−ζ) dv  +  g(ζ)·log((upper−ζ)/(lower−ζ))  [+ 2πi·g(ζ)]

The first integrand has a *removable* singularity at `v=ζ` (g analytic), so plain
adaptive quadrature handles it; the single complex `log` of the ratio carries the
branch cut (same invariant as `cell_hilbert`). For `Im ζ → 0⁺` the log limit
supplies the `+iπ g(ζ)` Plemelj term automatically; the explicit `2πi·g(ζ)` is
the Landau continuation added only when `ζ` drops below the real axis with
`Re ζ` inside the support (the growing-mode sheet). General-`g` sibling of the
`Z` overload: `hilbert(v->exp(-v^2)/√π, ζ; lower=-Inf, upper=Inf) == Z(ζ)`.
"""
function hilbert(g, ζ; lower, upper, rtol = 1.0e-9)
    gζ = g(ζ)
    reg = QuadGK.quadgk(v -> (g(v) - gζ) / (v - ζ), lower, upper; rtol)[1]
    return reg + gζ * _landau_logfac(ζ, lower, upper)
end

# Plemelj branch-cut-safe log ratio `log((hi−ζ)/(lo−ζ))` of the pole `1/(p−ζ)`
# integrated over `[lo,hi]`, plus the `2πi` Landau continuation when `ζ` sits on
# the growing-mode sheet (`Im ζ<0`, `Re ζ` in range).
@inline function _landau_logfac(ζ, lo, hi)
    logfac = log((hi - ζ) / (lo - ζ))
    return (imag(ζ) < 0 && lo < real(ζ) < hi) ? logfac + 2π * im : logfac
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
