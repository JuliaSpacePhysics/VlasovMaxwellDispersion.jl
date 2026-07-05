# The exact per-cell closed form (`hilbert_pwpoly`×`perp_pwpoly` on `fit.coeffs`) is the documented optimization
# over the adaptive-quadrature coupled path; not yet wired.

"""
    GridVDF(vperp, vpar, f; method=NonnegBSpline{3}(), regime=NonRelativistic())

Tabulated gyrotropic VDF: `f[i,j] = f₀(vperp[i], vpar[j])` on an ascending grid (`vperp[1] ≥ 0`).
A `method` turns the grid into a **tensor** `f₀` (knots + per-cell power coeffs) — a complex-evaluable `f₀` with analytic `∂∥/∂⊥`.

Fit methods (`fit_grid`): `NonnegBSpline` (default; positivity-preserving two-pass NNLS B-spline with `f ≥ 0`) or
`BicubicHermite` (local C¹ interpolation, O(N), no positivity guard).
The fit is renormalized to `∫d³p f₀ = 1`.
"""
struct GridVDF{V, C, F <: TensorSplineFit} <: AbstractVDF
    vpar::V
    vperp::V
    fit::F
    coupled::C          # CoupledVDF wrapping the complex-evaluable spline
end

regime(d::GridVDF) = regime(d.coupled)

function GridVDF(vperp, vpara, f; rtol = 1.0e-3, method = nothing, regime = NonRelativistic())
    method = @something(method, NonnegBSpline{3}(; rtol, maxknots_para = length(vpara), maxknots_perp = length(vperp)))
    # rescale as a tiny-valued grid (e.g. exp(-μγ)~1e-18) would otherwise underflow the fit to all-zeros
    scale = maximum(abs, f)
    iszero(scale) && throw(ArgumentError("GridVDF: f is all zeros"))
    fit = fit_grid(method, vperp, vpara, f ./ scale)   # (method, vperp, vpar, F[perp,par])
    fit.coeffs ./= _fit_d3p(fit)
    dgrad = (v, u) -> _grad2(fit, v, u)
    para = promote(float(fit.knots_para[1]), float(fit.knots_para[end]))
    perp = oftype(para[2], fit.knots_perp[1]), oftype(para[2], fit.knots_perp[end])
    cpl = CoupledVDF(fit, dgrad, para, perp, one(para[1]), regime)
    return GridVDF(vpara, vperp, fit, cpl)
end

# A tabulated f₀ is ZERO outside its sampled support — never the bicubic's cubic
# extrapolation (which diverges and breaks the relativistic Plemelj, where ∇f₀ is
# probed at the complex pole ζ far off-grid). Cell chosen by Re; outside ⇒ 0.
@inline function _insupport(fit::TensorSplineFit, v, u)
    return fit.knots_perp[1] <= real(v) <= fit.knots_perp[end] &&
        fit.knots_para[1] <= real(u) <= fit.knots_para[end]
end

# ∫d³p f₀ = 2π∫∫ p⊥ f₀ dp⊥ dv∥ in closed form from the cell coeffs (cylindrical
# weight p⊥ = knots_perp[i]+s⊥). cell[A,B] is s⊥^{A-1} s∥^{B-1}: A is the perp
# (first) axis, B the parallel.
# Used to normalize the fit
function _fit_d3p(fit::TensorSplineFit)
    kq, kp, c = fit.knots_perp, fit.knots_para, fit.coeffs
    acc = 0.0
    @inbounds for i in 1:(length(kq) - 1), j in 1:(length(kp) - 1)
        hq, hp, wl = kq[i + 1] - kq[i], kp[j + 1] - kp[j], kq[i]
        cell = c[i, j]
        for A in axes(cell, 1), B in axes(cell, 2)
            acc += cell[A, B] * (hp^B / B) * (wl * hq^A / A + hq^(A + 1) / (A + 1))
        end
    end
    return 2π * acc
end

# Exact parallel moments make the GridVDF far faster than the generic
# coupled path: f₀ IS piecewise-polynomial, so the inner Landau H∥ closes per cell
# (`cell_hilbert_landau`) instead of adaptive QuadGK.
# Relativistic grids route through the coupled (p⊥,p∥) path: it integrates exactly the
# grid rectangle (only p∥ goes complex, at the near-axis poles), so the spline is never
# probed off-grid.
function contribution(d::GridVDF, s, ω, k; closure::IntegralClosure = HarmonicSum())
    if closure isa HarmonicSum && regime(d) isa NonRelativistic
        return _grid_contribution(d, s, complex(float(ω)), k)
    end
    return contribution(d.coupled, s, ω, k; closure)
end

function _grid_contribution(d::GridVDF, s, ω, k; rtol = 1.0e-6)
    fit = d.fit
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    a = kperp / Ω
    # nmax from the perp scale (mean p⊥² over the fitted), as in CoupledVDF.
    p⊥²_mean = 2π * QuadGK.quadgk(
        v -> v^3 * _fit_par_integral(fit, v), zero(fit.knots_perp[end]), fit.knots_perp[end]; rtol = 1.0e-7
    )[1]
    nmax = nmax_bessel(a^2 * abs(p⊥²_mean) / 2)
    f = n -> _grid_harmonic(n, fit, ω, Ω, kz, a)
    χ = converge(f; nmax, rtol)
    return (s.Pi2 / ω^2) * _antisymmat(χ)
end

# ∫ f₀(v,u) du over the parallel support at fixed perp v (for the perp-scale estimate).
# cell[A,B] is s⊥^{A-1} s∥^{B-1}: t=s⊥ fixes the perp cell i; integrate s∥ (B axis).
@inline function _fit_par_integral(fit::TensorSplineFit, v)
    i = _cell(fit.knots_perp, v)
    t = v - fit.knots_perp[i]
    kp = fit.knots_para
    acc = 0.0
    @inbounds for j in 1:(length(kp) - 1)
        h = kp[j + 1] - kp[j]
        cell = fit.coeffs[i, j]
        for A in axes(cell, 1), B in axes(cell, 2)
            acc += cell[A, B] * t^(A - 1) * h^B / B
        end
    end
    return acc
end

# Local poly with coeffs `a` (ascending, in s=u-vl) → absolute-u monomial coeffs
# (same length): bₘ = Σ_{k≥m} a[k] C(k,m) (-vl)^{k-m}.
@inline function _shift_to_abs(a::SVector{L}, vl) where {L}
    return SVector{L}(
        ntuple(L) do m1
            m = m1 - 1
            s = zero(eltype(a))
            @inbounds for k in m:(L - 1)
                s += a[k + 1] * binomial(k, m) * (-vl)^(k - m)
            end
            s
        end
    )
end

# `cell_hilbert_landau` with the per-cell `log((vr-ζ)/(vl-ζ))` and Landau flag
# precomputed (both depend only on the cell and ζ, not on the coeffs) — and it
# reuses the Horner remainder `pζ=p(ζ)` for the Landau term instead of a second
# `_polyval`. Hot-loop inner kernel for the parallel moments.
@inline function _cellH(coeffs, vl, vr, ζ, logr, landau)
    m = length(coeffs)
    T = complex(promote_type(eltype(coeffs), typeof(ζ)))
    pζ = convert(T, coeffs[m])
    poly = zero(T)
    @inbounds for k in (m - 1):-1:1
        d = k - 1
        poly += pζ * (vr^(d + 1) - vl^(d + 1)) / (d + 1)
        pζ = coeffs[k] + ζ * pζ
    end
    h = poly + pζ * logr
    return landau ? h + 2π * im * pζ : h
end

# Per-perp-cell t-POLYNOMIAL coefficients of the 5 parallel moments at pole ζ.
# Key identity: the p⊥-slice coeffs are polynomials in t (= p⊥ − knots_perp[i]),
# and `cell_hilbert` is linear in them ⇒ each moment z(t) is a polynomial in t —
# F moments (∂⊥ slice) deg NP-2, T moments (∂∥ slice) deg NP-1. We compute it ONCE
# per perp cell (and the perp-para-cell log ONCE per harmonic), then the p⊥ quadrature
# only evaluates that polynomial per node instead of re-summing `cell_hilbert` over
# every parallel cell. Exact; −1/kz folds the resonance kz.
# cell[A,B] is s⊥^{A-1} s∥^{B-1}: A is the perp (t) axis, B the parallel (Hilbert) axis.
@inline function _grid_parmoment_polys(fit::TensorSplineFit, i, ζ)
    kp = fit.knots_para
    c = fit.coeffs
    coeff_Type = eltype(c)
    T = complex(promote_type(eltype(coeff_Type), typeof(ζ)))
    NP, NQ = size(coeff_Type)
    MF0 = MF1 = MF2 = zero(SVector{NP - 1, T})   # t^0..t^{NP-2}
    MT0 = MT1 = zero(SVector{NP, T})             # t^0..t^{NP-1}
    @inbounds for j in 1:(length(kp) - 1)
        vl, vr = kp[j], kp[j + 1]
        cell = c[i, j]
        logr = log((vr - ζ) / (vl - ζ))
        landau = imag(ζ) < 0 && _pole_in_cell(vl, vr, ζ)
        # ∂⊥ slice: t-power b ⇐ row b+2, weight (b+1); s∥-poly P (length NQ).
        for b in 0:(NP - 2)
            w = b + 1
            row = b + 2
            P = _shift_to_abs(SVector(ntuple(B -> w * cell[row, B], Val(NQ))), vl)
            hF0 = _cellH(P, vl, vr, ζ, logr, landau)
            hF1 = _cellH(vcat(SVector(zero(eltype(P))), P), vl, vr, ζ, logr, landau)
            hF2 = _cellH(vcat(SVector(zero(eltype(P)), zero(eltype(P))), P), vl, vr, ζ, logr, landau)
            MF0 = setindex(MF0, MF0[b + 1] + hF0, b + 1)
            MF1 = setindex(MF1, MF1[b + 1] + hF1, b + 1)
            MF2 = setindex(MF2, MF2[b + 1] + hF2, b + 1)
        end
        # ∂∥ slice: t-power b ⇐ row b+1; s∥-deriv poly Q (length NQ-1, coeff
        # of s∥^{m-1} is m·cell[row,m+1]).
        for b in 0:(NP - 1)
            row = b + 1
            Q = _shift_to_abs(SVector(ntuple(m -> m * cell[row, m + 1], Val(NQ - 1))), vl)
            hT0 = _cellH(Q, vl, vr, ζ, logr, landau)
            hT1 = _cellH(vcat(SVector(zero(eltype(Q))), Q), vl, vr, ζ, logr, landau)
            MT0 = setindex(MT0, MT0[b + 1] + hT0, b + 1)
            MT1 = setindex(MT1, MT1[b + 1] + hT1, b + 1)
        end
    end
    return (MF0, MF1, MF2, MT0, MT1)
end

# One cyclotron harmonic: loop perp cells, precompute the parallel-moment
# t-polynomials once per cell, then a smooth p⊥ QuadGK whose integrand only
# evaluates those + Bessel weights.
function _grid_harmonic(n, fit::TensorSplineFit, ω, Ω, kz, a)
    ζ = (ω - n * Ω) / kz
    kq = fit.knots_perp
    acc = zero(SVector{6, ComplexF64})
    for i in 1:(length(kq) - 1)
        wl, wr = kq[i], kq[i + 1]
        MF0c, MF1c, MF2c, MT0c, MT1c = _grid_parmoment_polys(fit, i, ζ)
        integ = v -> begin
            t = v - wl
            M = (
                evalpoly(t, MF0c.data), evalpoly(t, MF1c.data), evalpoly(t, MF2c.data),
                evalpoly(t, MT0c.data), evalpoly(t, MT1c.data),
            )
            _In_block(M, (-1 / kz), _perp_Bessel_bilinear(n, a, v), v, ω, kz, n * Ω)
        end
        # The Bessel weight J_n(a v) has v-wavelength ≈ π/a; adaptive QuadGK over a
        # cell spanning many wavelengths can't resolve it. Pre-split the cell at the
        # oscillation scale (≥2 panels/wavelength) and let QuadGK adapt within each.
        seg = first(_quadgk_osc(integ, wl, wr, a))
        acc = acc .+ seg
    end
    return acc
end

# Adaptive QuadGK with oscillation-scale breakpoints for the Bessel kernel.
@inline function _quadgk_osc(f, wl, wr, a)
    nb = ceil(Int, abs(a) * (wr - wl) / (π / 2))
    pts = nb <= 1 ? (wl, wr) : range(wl, wr; length = nb + 1)
    return QuadGK.quadgk(f, pts...; rtol = 1.0e-6, norm = NORM)
end
