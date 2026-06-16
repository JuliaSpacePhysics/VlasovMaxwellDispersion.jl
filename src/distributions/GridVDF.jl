# The exact per-cell closed form (`hilbert_pwpoly`×`perp_pwpoly` on `fit.coeffs`) is the documented optimization
# over the adaptive-quadrature coupled path; not yet wired.

"""
    GridVDF(vpar, vperp, f; method=NonnegBSpline{3}())

Tabulated gyrotropic VDF: `f[i,j] = f₀(vpar[i], vperp[j])` on an ascending grid
(`vperp[1] ≥ 0`). A `method` turns the grid into a **tensor** `f₀` (knots + per-cell power coeffs) — 
a complex-evaluable `f₀` with analytic `∂∥/∂⊥`.

Fit methods (`fit_grid`): `NonnegBSpline` (default; positivity-preserving two-pass NNLS B-spline with `f ≥ 0`) or
`BicubicHermite` (local C¹ interpolation, O(N), no positivity guard).
The fit is renormalized to `∫d³p f₀ = 1`.
"""
struct GridVDF{C,F<:TensorSplineFit} <: AbstractVDF
    vpar::Vector{Float64}
    vperp::Vector{Float64}
    fit::F
    coupled::C          # CoupledVDF wrapping the complex-evaluable spline
end

function GridVDF(vpar, vperp, f; tol=1.0e-3, method=nothing)
    method = @something(method, NonnegBSpline{3}(; tol, maxknots_par=length(vpar), maxknots_perp=length(vperp)))
    vp, vq = collect(float.(vpar)), collect(float.(vperp))
    # rescale as a tiny-valued grid (e.g. exp(-μγ)~1e-18) would otherwise underflow the fit to all-zeros
    scale = maximum(abs, f)
    iszero(scale) && throw(ArgumentError("GridVDF: f is all zeros"))
    fit = fit_grid(method, vp, vq, f ./ scale)
    fit.coeffs ./= _fit_d3p(fit)
    dpar = (u, v) -> _fit_dpar(fit, u, v)
    dperp = (u, v) -> _fit_dperp(fit, u, v)
    cpl = CoupledVDF(
        fit; parlower=fit.knots_par[1], parupper=fit.knots_par[end],
        perpupper=fit.knots_perp[end], dpar, dperp, normalize=false
    )
    GridVDF(vp, vq, fit, cpl)
end

# A tabulated f₀ is ZERO outside its sampled support — never the bicubic's cubic
# extrapolation (which diverges and breaks the relativistic Plemelj, where ∇f₀ is
# probed at the complex pole ζ far off-grid). Cell chosen by Re; outside ⇒ 0.
@inline function _insupport(fit::TensorSplineFit, u, v)
    fit.knots_par[1] <= real(u) <= fit.knots_par[end] &&
        fit.knots_perp[1] <= real(v) <= fit.knots_perp[end]
end


@inline function _fit_dpar(fit::TensorSplineFit, u, v)
    _insupport(fit, u, v) || return zero(complex(promote_type(typeof(float(u)), typeof(float(v)))))
    i, j = _cell(fit.knots_par, u), _cell(fit.knots_perp, v)
    _dpolyval2_s(fit.coeffs[i, j], u - fit.knots_par[i], v - fit.knots_perp[j])
end
@inline function _fit_dperp(fit::TensorSplineFit, u, v)
    _insupport(fit, u, v) || return zero(complex(promote_type(typeof(float(u)), typeof(float(v)))))
    i, j = _cell(fit.knots_par, u), _cell(fit.knots_perp, v)
    _dpolyval2_t(fit.coeffs[i, j], u - fit.knots_par[i], v - fit.knots_perp[j])
end

# ∫d³p f₀ = 2π∫∫ p⊥ f₀ dv∥dp⊥ in closed form from the cell coeffs (cylindrical
# weight p⊥ = knots_perp[j]+t). Used to normalize the fit so every VDF shares the
# ∫d³p=1 convention — `fit_grid` itself only preserves shape.
function _fit_d3p(fit::TensorSplineFit)
    kp, kq, c = fit.knots_par, fit.knots_perp, fit.coeffs
    acc = 0.0
    @inbounds for i in 1:(length(kp)-1), j in 1:(length(kq)-1)
        hp, hq, wl = kp[i+1] - kp[i], kq[j+1] - kq[j], kq[j]
        cell = c[i, j]
        for A in axes(cell, 1), B in axes(cell, 2)
            acc += cell[A, B] * (hp^A / A) * (wl * hq^B / B + hq^(B + 1) / (B + 1))
        end
    end
    2π * acc
end

# Exact parallel moments make the GridVDF far faster than the generic
# coupled path: f₀ IS piecewise-polynomial, so the inner Landau H∥ closes per cell
# (`cell_hilbert_landau`) instead of adaptive QuadGK.
function contribution(d::GridVDF, s::Species, ω, k; closure::IntegralClosure=HarmonicSum())
    iszero(perp(k)) &&
        throw(ArgumentError("GridVDF: magnetized EM tensor needs kperp≠0 (oblique)"))
    if closure isa HarmonicSum
        Regime(s) isa NonRelativistic && return _grid_contribution(d, s, complex(float(ω)), k)
        Regime(s) isa Relativistic && return _grid_contribution_rel(d, s, complex(float(ω)), k)
    end
    return contribution(d.coupled, s, ω, k; closure)
end

# Relativistic grid path. The grid tabulates f₀ on a (p∥,p⊥) RECTANGLE, but the
# (γ,p∥) integral sweeps constant-γ half-circles p∥²+p⊥²=γ²−1. Feeding the
# straight to the coupled-rel evaluator fails (NaN): its γmax=√(1+p∥max²+p⊥max²)
# reaches the rectangle's CORNER, where the half-circle pokes outside the grid
# (p⊥>perphi) and the extrapolates to garbage. Fix: integrate only up to
# the largest γ whose whole half-circle fits inside the grid — radius R=min(p∥
# half-width, perphi) ⇒ γmax=√(1+R²). f₀ decays, so the clipped corner is
# negligible. Otherwise this reuses the validated `_coupled_harmonic_rel` (fixed
# GL in (γ,p∥), inner single-pole Plemelj) with the analytic ∇f₀.
function _grid_contribution_rel(d::GridVDF, s::Species, ω, k)
    cpl = d.coupled
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    a = kperp / Ω
    R = min(min(cpl.parhi, -cpl.parlo), cpl.perphi)   # half-circle radius that fits the grid
    γmax = sqrt(1 + R^2)
    nmax = nmax_bessel(a^2 * R^2 / 2)
    f = n -> _grid_harmonic_rel(n, cpl, ω, Ω, kz, a, γmax)
    χ = converge(f, 1, 1.0e-6; nmax)
    χ = χ .+ _ee33(_bernstein_rel(cpl, ω, γmax))   # non-resonant Bernstein (drops from the n-sum)
    return SMatrix{3,3,ComplexF64}((s.Pi2 / ω^2) * χ)
end

# Relativistic harmonic. Like `_coupled_harmonic_rel` but the
# Plemelj pole subtraction — which evaluates ∇f₀ at p∥=ζ, hence p⊥=√(γ²−1−ζ²) —
# is done ONLY when the resonance is in range (|Re ζ|<umax ⇒ p⊥ real,). 
# When ζ is off the real momentum range the integrand g/(p∥−ζ) is already
# smooth, so we integrate it directly and never probe the spine off-grid (where
# p⊥ is imaginary and the local polynomial has no analytic continuation).
function _grid_harmonic_rel(n, d, ω, Ω, kz, a, γmax)
    gn, gw = _GLγ
    pn, pw = _GLp
    acc = zero(SMatrix{3,3,ComplexF64})
    for ig in eachindex(gn)
        γ, wγ = _rescale(gn[ig], gw[ig], one(real(ω)), γmax)
        umax = sqrt(γ^2 - 1)
        ζ = (γ * ω - n * Ω) / kz
        if -umax < real(ζ) < umax
            nζ = _rel_integrand(ζ, γ, n, a, ω, kz, d)   # p⊥(ζ) real ⇒ bicubic valid
            reg = zero(nζ)
            for ip in eachindex(pn)
                u, wu = _rescale(pn[ip], pw[ip], -umax, umax)
                reg = reg .+ wu .* ((_rel_integrand(u, γ, n, a, ω, kz, d) .- nζ) ./ (u - ζ))
            end
            inner = reg .+ nζ .* _landau_logfac(ζ, -umax, umax)
        else
            inner = zero(SMatrix{3,3,ComplexF64})       # smooth: integrate directly
            for ip in eachindex(pn)
                u, wu = _rescale(pn[ip], pw[ip], -umax, umax)
                inner = inner .+ wu .* (_rel_integrand(u, γ, n, a, ω, kz, d) ./ (u - ζ))
            end
        end
        acc = acc .+ wγ .* ((-1 / kz) .* inner)
    end
    return acc
end

function _grid_contribution(d::GridVDF, s::Species, ω, k)
    fit = d.fit
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    a = kperp / Ω
    # nmax from the perp scale (mean p⊥² over the fitted), as in CoupledVDF.
    p⊥²_mean = 2π * QuadGK.quadgk(
        v -> v^3 * _fit_par_integral(fit, v), zero(fit.knots_perp[end]), fit.knots_perp[end]; rtol=1.0e-7
    )[1]
    nmax = nmax_bessel(a^2 * abs(p⊥²_mean) / 2)
    f = n -> _grid_harmonic(n, fit, ω, Ω, kz, kperp, a)
    χ = converge(f, 1, 1.0e-6; nmax)
    return SMatrix{3,3,ComplexF64}((s.Pi2 / ω^2) * χ)
end

# ∫ f₀(u,v) du over the parallel support at fixed v (for the perp-scale estimate).
@inline function _fit_par_integral(fit::TensorSplineFit, v)
    j = _cell(fit.knots_perp, v)
    t = v - fit.knots_perp[j]
    kp = fit.knots_par
    acc = 0.0
    @inbounds for i in 1:(length(kp)-1)
        h = kp[i+1] - kp[i]
        cell = fit.coeffs[i, j]
        for A in axes(cell, 1), B in axes(cell, 2)
            acc += cell[A, B] * t^(B - 1) * h^A / A
        end
    end
    acc
end

# Local poly with coeffs `a` (ascending, in s=u-vl) → absolute-u monomial coeffs
# (same length): bₘ = Σ_{k≥m} a[k] C(k,m) (-vl)^{k-m}.
@inline function _shift_to_abs(a::SVector{L}, vl) where {L}
    SVector{L}(ntuple(L) do m1
        m = m1 - 1
        s = zero(eltype(a))
        @inbounds for k in m:(L-1)
            s += a[k+1] * binomial(k, m) * (-vl)^(k - m)
        end
        s
    end)
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
    @inbounds for k in (m-1):-1:1
        d = k - 1
        poly += pζ * (vr^(d + 1) - vl^(d + 1)) / (d + 1)
        pζ = coeffs[k] + ζ * pζ
    end
    h = poly + pζ * logr
    return landau ? h + 2π * im * pζ : h
end

# Per-perp-cell t-POLYNOMIAL coefficients of the 5 parallel moments at pole ζ.
# Key identity: the p⊥-slice coeffs are polynomials in t (= p⊥ − knots_perp[j]),
# and `cell_hilbert` is linear in them ⇒ each moment z(t) is a polynomial in t —
# F moments (∂⊥ slice) deg NQ-2, T moments (∂∥ slice) deg NQ-1. We compute it ONCE
# per perp cell (and the per-par-cell log ONCE per harmonic), then the p⊥ quadrature
# only evaluates that polynomial per node instead of re-summing `cell_hilbert` over
# every parallel cell. Exact; −1/kz folds the resonance kz.
@inline function _grid_parmoment_polys(fit::TensorSplineFit{Tc,NP,NQ}, j, ζ, kz) where {Tc,NP,NQ}
    kp = fit.knots_par
    c = fit.coeffs
    T = complex(promote_type(Tc, typeof(ζ)))
    z0F = z1F = z2F = zero(SVector{NQ - 1,T})   # t^0..t^{NQ-2}
    z0T = z1T = zero(SVector{NQ,T})             # t^0..t^{NQ-1}
    @inbounds for i in 1:(length(kp)-1)
        vl, vr = kp[i], kp[i+1]
        cell = c[i, j]
        logr = log((vr - ζ) / (vl - ζ))
        landau = imag(ζ) < 0 && _pole_in_cell(vl, vr, ζ)
        # ∂⊥ slice: t-power b ⇐ column b+2, weight (b+1); s∥-poly P (length NP).
        for b in 0:(NQ-2)
            w = b + 1
            col = b + 2
            P = _shift_to_abs(SVector(ntuple(A -> w * cell[A, col], Val(NP))), vl)
            hF0 = _cellH(P, vl, vr, ζ, logr, landau)
            hF1 = _cellH(vcat(SVector(zero(eltype(P))), P), vl, vr, ζ, logr, landau)
            hF2 = _cellH(vcat(SVector(zero(eltype(P)), zero(eltype(P))), P), vl, vr, ζ, logr, landau)
            z0F = setindex(z0F, z0F[b+1] + hF0, b + 1)
            z1F = setindex(z1F, z1F[b+1] + hF1, b + 1)
            z2F = setindex(z2F, z2F[b+1] + hF2, b + 1)
        end
        # ∂∥ slice: t-power b ⇐ column b+1; s∥-deriv poly Q (length NP-1, coeff
        # of s^{m-1} is m·cell[m+1,col]).
        for b in 0:(NQ-1)
            col = b + 1
            Q = _shift_to_abs(SVector(ntuple(m -> m * cell[m+1, col], Val(NP - 1))), vl)
            hT0 = _cellH(Q, vl, vr, ζ, logr, landau)
            hT1 = _cellH(vcat(SVector(zero(eltype(Q))), Q), vl, vr, ζ, logr, landau)
            z0T = setindex(z0T, z0T[b+1] + hT0, b + 1)
            z1T = setindex(z1T, z1T[b+1] + hT1, b + 1)
        end
    end
    f = -1 / kz
    return (f * z0F, f * z1F, f * z2F, f * z0T, f * z1T)
end

# One cyclotron harmonic: loop perp cells, precompute the parallel-moment
# t-polynomials once per cell, then a smooth p⊥ QuadGK whose integrand only
# evaluates those + Bessel weights.
function _grid_harmonic(n, fit::TensorSplineFit, ω, Ω, kz, kperp, a)
    ζ = (ω - n * Ω) / kz
    kq = fit.knots_perp
    acc = zero(SMatrix{3,3,ComplexF64})
    for j in 1:(length(kq)-1)
        wl, wr = kq[j], kq[j+1]
        z0Fc, z1Fc, z2Fc, z0Tc, z1Tc = _grid_parmoment_polys(fit, j, ζ, kz)
        seg = first(QuadGK.quadgk(wl, wr; rtol=1.0e-6, norm=x -> maximum(abs, x)) do v
            t = v - wl
            z = (evalpoly(t, z0Fc.data), evalpoly(t, z1Fc.data), evalpoly(t, z2Fc.data),
                evalpoly(t, z0Tc.data), evalpoly(t, z1Tc.data))
            return _chi_mblock(z, _perp_bessel_moments(n, a, v), ω, kz, kperp, n / a)
        end)
        acc = acc .+ seg
    end
    return acc
end
