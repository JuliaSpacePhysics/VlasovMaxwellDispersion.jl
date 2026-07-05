# Standalone Landau-continued Cauchy transform evaluator (sinc/cot method).
# Self-contained — no dependencies beyond Base. See README.md for the derivation.
#
# Composable interface:  cauchy_landau(g, ζ[, alg])
#   alg = SincCot(map = SinhMap(uc, S), h = 0.2, window = 7.0)
# The cot pole correction is MAP-INDEPENDENT — for any analytic bijection
# ψ: ℝ → ℝ the residue of G(t)/(ψ(t)−ζ) at t* = ψ⁻¹(ζ) is exactly g(ζ)
# (Jacobian cancels) — so maps compose freely with the same correction. A map
# only changes where nodes land (tail compression, feature centering) and the
# strip geometry that sets the error rate e^(−2πd_t/h).

abstract type AbstractMap end

"""ψ(t) = uc + S·sinh(t): tails compressed logarithmically — the default.
Handles algebraic (kappa-like) decay with a logarithmic node window."""
struct SinhMap{T} <: AbstractMap
    uc::T
    S::T
end
SinhMap(; uc = 0.0, S = 1.0) = SinhMap(promote(uc, S)...)
mapto(m::SinhMap, t) = m.uc + m.S * sinh(t)
mapjac(m::SinhMap, t) = m.S * cosh(t)
mapinv(m::SinhMap, ζ) = asinh((ζ - m.uc) / m.S)

"""ψ(t) = uc + S·t: plain trapezoid on ℝ. Optimal when g is entire with fast
(e.g. Gaussian) decay — the window must cover the support in units of S, so
heavy tails want [`SinhMap`](@ref) instead."""
struct LinearMap{T} <: AbstractMap
    uc::T
    S::T
end
LinearMap(; uc = 0.0, S = 1.0) = LinearMap(promote(uc, S)...)
mapto(m::LinearMap, t) = m.uc + m.S * t
mapjac(m::LinearMap, t) = m.S
mapinv(m::LinearMap, ζ) = (ζ - m.uc) / m.S

"""
    SincCot(; map=SinhMap(), h=0.2, window=7.0)

Mapped-trapezoid + cotangent pole correction:

    C⁺[g](ζ) = h·Σⱼ G(tⱼ)/(ψ(tⱼ)−ζ) + π·g(ζ)·(cot(π·t*/h) + i),
    G = g(ψ)·ψ′,  t* = ψ⁻¹(ζ),  tⱼ = j·h,  |j| ≤ window/h.

Error ~ e^(−2πd/h), d = strip half-width of G in t (h=0.2 → ~1e-9 for
thermal-like g). Nodes are exact multiples of h — required by the cot phase.
"""
struct SincCot{M <: AbstractMap, T <: Real}
    map::M
    h::T
    window::T
end
SincCot(; map = SinhMap(), h = 0.2, window = 7.0) = SincCot(map, promote(h, window)...)

"""
    cauchy_landau(g, ζ,  alg = SincCot())
    cauchy_landau(g, ζs, alg = SincCot())

Landau-continued Cauchy transform of a strip-analytic decaying `g`,

    C⁺[g](ζ) = ∫_ℝ g(u)/(u−ζ) du   continued analytically from Im ζ > 0
             (Im ζ < 0: integral + 2πi·g(ζ);  Im ζ = 0: PV + iπ·g(ζ)).

The vector method shares the g samples across all `ζs` — the marginal cost of
a ladder point is one `g(ζ)` and one `cot`. `g` may return any value with
linear arithmetic (scalars, static vectors); it must not overflow at large
real |u|.
"""
cauchy_landau(g, ζ::Number, alg::SincCot = SincCot()) = only(cauchy_landau(g, (ζ,), alg))

function cauchy_landau(g, ζs, alg::SincCot = SincCot())
    (; map, h, window) = alg
    jmax = floor(Int, window / h)
    u0 = mapto(map, zero(h))
    acc = [zero(h * g(u0) / (u0 - ζ)) for ζ in ζs]
    for j in (-jmax):jmax
        t = j * h                          # nodes MUST be exact multiples of h:
        u = mapto(map, t)                  # an offset grid breaks the cot phase
        w = h * mapjac(map, t)
        gu = g(u)
        for (i, ζ) in pairs(ζs)
            acc[i] += w * gu / (u - ζ)
        end
    end
    return [acc[i] + π * g(ζ) * _cot_i(π * mapinv(map, ζ) / h) for (i, ζ) in pairs(ζs)]
end

# cot(w) + i, saturated: cot → ∓i as Im w → ±∞ but overflows past |Im w| ≈ 700;
# beyond 20 the residual is < 4e-18, so return the limit (0 or 2i) exactly.
_cot_i(w) = abs(imag(w)) > 20 ? complex(zero(real(w)), imag(w) > 0 ? zero(real(w)) : 2 * one(real(w))) : cot(w) + im
