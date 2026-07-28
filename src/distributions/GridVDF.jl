"""
    GridVDF(vperp, vpar, f; method=NonnegBSpline{3}(), regime=NonRelativistic())

Tabulated gyrotropic VDF: `f[i,j] = f₀(vperp[i], vpar[j])` on an ascending grid (`vperp[1] ≥ 0`).
A `method` projects the grid onto a compactly supported, complex-evaluable `f₀`.

Built-in methods return tensor piecewise polynomials: `NonnegBSpline` (default;
positivity-preserving two-pass NNLS B-spline with `f ≥ 0`) or `BicubicHermite`
(local C¹ interpolation, O(N), no positivity guard).
Susceptibility is normalized by the projected density.
"""
struct GridVDF{F,C} <: AbstractVDF
    fit::F
    coupled::C
end

regime(d::GridVDF) = regime(d.coupled)

_vgrid(v) = eltype(v) <: Real ? v : map(velocity, v) # Unitful support

function GridVDF(vperp, vpara, f; rtol = 1.0e-3, method = nothing, regime = NonRelativistic())
    vperp, vpara = _vgrid(vperp), _vgrid(vpara)
    method = @something(method, NonnegBSpline{3}(; rtol, maxknots_para = length(vpara), maxknots_perp = length(vperp)))
    # rescale as a tiny-valued grid (e.g. exp(-μγ)~1e-18) would otherwise underflow the fit to all-zeros
    scale = maximum(abs, f)
    iszero(scale) && throw(ArgumentError("GridVDF: f is all zeros"))
    fit = fit_grid(method, vperp, vpara, f ./ scale)
    n = 2π * moment(fit, (1, 0))
    iszero(n) && throw(ArgumentError("GridVDF: projection has zero density"))
    dgrad = (v, u) -> _grad2(fit, v, u)
    lims = domain(fit)
    para = promote(float(lims.para[1]), float(lims.para[2]))
    perp = oftype(para[2], lims.perp[1]), oftype(para[2], lims.perp[2])
    cpl = CoupledVDF(erase_f2(fit, para[2]), erase_g2(dgrad, para[2]), nothing, para, perp, regime)
    cache = if regime isa NonRelativistic
        p2m = 2π * moment(fit, (3, 0)) / n
        (; n = oftype(para[1], n), pperp2_mean = oftype(para[1], p2m))
    else
        (; n = oftype(para[1], n), bernstein33 = _bernstein_rel(cpl))
    end
    return GridVDF(fit, PreparedVDF(cpl, cache))
end

# Exact parallel moments make the tensor projection faster than the generic coupled path.
# Relativistic grids route through the coupled (p⊥,p∥) path: it integrates exactly the
# grid rectangle (only p∥ goes complex, at the near-axis poles), so the spline is never
# probed off-grid.
function contribution(d::GridVDF, s, ω, k; closure::IntegralClosure = HarmonicSum())
    if closure isa HarmonicSum && regime(d) isa NonRelativistic
        return _grid_contribution(d.fit, d, s, complex(float(ω)), k)
    end
    return contribution(d.coupled, s, ω, k; closure)
end

_grid_contribution(fit, d, s, ω, k; kw...) =
    contribution(d.coupled, s, ω, k; closure = HarmonicSum(), kw...)

function _grid_contribution(fit::TensorSplineFit, d, s, ω, k; rtol = 1.0e-6)
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    a = kperp / Ω
    nmax = nmax_bessel(a^2 * abs(d.coupled.cache.pperp2_mean) / 2)
    f = n -> _grid_harmonic(n, fit, ω, Ω, kz, a)
    χ = converge(f; nmax, rtol)
    return (s.Pi2 / (ω^2 * d.coupled.cache.n)) * _antisymmat(χ)
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

# Cell Hilbert integral with log and Landau jump shared across all moments.
@inline function _cell_hilbert(coeffs, vl, vr, ζ, logr, jump)
    m = length(coeffs)
    T = complex(promote_type(eltype(coeffs), typeof(ζ)))
    pζ = convert(T, coeffs[m])
    poly = zero(T)
    @inbounds for k in (m - 1):-1:1
        d = k - 1
        poly += pζ * (vr^(d + 1) - vl^(d + 1)) / (d + 1)
        pζ = coeffs[k] + ζ * pζ
    end
    return poly + pζ * (logr + jump)
end

# Per-perp-cell polynomial coefficients of the five parallel moments at pole ζ.
# Key identity: the p⊥-slice coeffs are polynomials in t (= p⊥ − knots_perp[i]),
# and `_cell_hilbert` is linear in them ⇒ each moment z(t) is a polynomial in t —
# F moments (∂⊥ slice) deg NP-2, T moments (∂∥ slice) deg NP-1. We compute it ONCE
# per perp cell (and the perp-para-cell log ONCE per harmonic), then the p⊥ quadrature
# only evaluates that polynomial per node instead of re-summing `_cell_hilbert` over
# every parallel cell. Exact; −1/kz folds the resonance kz.
# cell[A,B] is s⊥^{A-1} s∥^{B-1}: A is the perp (t) axis, B the parallel (Hilbert) axis.
@inline function _grid_parmoment_polys(fit::TensorSplineFit, i, ζ, σ)
    kp = fit.knots_para
    c = fit.coeffs
    Cell = eltype(c)
    T = complex(promote_type(eltype(Cell), typeof(ζ)))
    NP, NQ = size(Cell)
    MF0 = MF1 = MF2 = zero(SVector{NP - 1, T})   # t^0..t^{NP-2}
    MT0 = MT1 = zero(SVector{NP, T})             # t^0..t^{NP-1}
    @inbounds for j in 1:(length(kp) - 1)
        vl, vr = kp[j], kp[j + 1]
        cell = c[i, j]
        logr = log((vr - ζ) / (vl - ζ))
        jump = σ * imag(ζ) < 0 && vl < real(ζ) <= vr ? σ * 2π * im : zero(σ * 2π * im)
        # ∂⊥ slice: t-power b ⇐ row b+2, weight (b+1); s∥-poly P (length NQ).
        for b in 0:(NP - 2)
            P = _shift_to_abs(SVector(ntuple(B -> (b + 1) * cell[b + 2, B], Val(NQ))), vl)
            hF0 = _cell_hilbert(P, vl, vr, ζ, logr, jump)
            hF1 = _cell_hilbert(vcat(SVector(zero(eltype(P))), P), vl, vr, ζ, logr, jump)
            hF2 = _cell_hilbert(vcat(SVector(zero(eltype(P)), zero(eltype(P))), P), vl, vr, ζ, logr, jump)
            MF0 = setindex(MF0, MF0[b + 1] + hF0, b + 1)
            MF1 = setindex(MF1, MF1[b + 1] + hF1, b + 1)
            MF2 = setindex(MF2, MF2[b + 1] + hF2, b + 1)
        end
        # ∂∥ slice: t-power b ⇐ row b+1; s∥-deriv poly Q (length NQ-1, coeff
        # of s∥^{m-1} is m·cell[row,m+1]).
        for b in 0:(NP - 1)
            Q = _shift_to_abs(SVector(ntuple(m -> m * cell[b + 1, m + 1], Val(NQ - 1))), vl)
            hT0 = _cell_hilbert(Q, vl, vr, ζ, logr, jump)
            hT1 = _cell_hilbert(vcat(SVector(zero(eltype(Q))), Q), vl, vr, ζ, logr, jump)
            MT0 = setindex(MT0, MT0[b + 1] + hT0, b + 1)
            MT1 = setindex(MT1, MT1[b + 1] + hT1, b + 1)
        end
    end
    return (MF0, MF1, MF2, MT0, MT1)
end

# kz=0 kernel: ∫p(u)du over the cell (coeffs ascending in absolute u).
@inline function _cell_moment(coeffs, vl, vr)
    acc = zero(promote_type(eltype(coeffs), typeof(vl)))
    @inbounds for k in eachindex(coeffs)
        acc += coeffs[k] * (vr^k - vl^k) / k
    end
    return acc
end

# kz=0 variant of `_grid_parmoment_polys`: the Hilbert kernel degenerates to plain
# moments (`_cell_moment`, ζ-free) — the 1/Δ_n weight is applied by the caller.
@inline function _grid_parmoment_polys0(fit::TensorSplineFit, i)
    kp = fit.knots_para
    c = fit.coeffs
    Cell = eltype(c)
    T = eltype(Cell)
    NP, NQ = size(Cell)
    MF0 = MF1 = MF2 = zero(SVector{NP - 1, T})
    MT0 = MT1 = zero(SVector{NP, T})
    z = zero(T)
    @inbounds for j in 1:(length(kp) - 1)
        vl, vr = kp[j], kp[j + 1]
        cell = c[i, j]
        for b in 0:(NP - 2)
            P = _shift_to_abs(SVector(ntuple(B -> (b + 1) * cell[b + 2, B], Val(NQ))), vl)
            MF0 = setindex(MF0, MF0[b + 1] + _cell_moment(P, vl, vr), b + 1)
            MF1 = setindex(MF1, MF1[b + 1] + _cell_moment(vcat(SVector(z), P), vl, vr), b + 1)
            MF2 = setindex(MF2, MF2[b + 1] + _cell_moment(vcat(SVector(z, z), P), vl, vr), b + 1)
        end
        for b in 0:(NP - 1)
            Q = _shift_to_abs(SVector(ntuple(m -> m * cell[b + 1, m + 1], Val(NQ - 1))), vl)
            MT0 = setindex(MT0, MT0[b + 1] + _cell_moment(Q, vl, vr), b + 1)
            MT1 = setindex(MT1, MT1[b + 1] + _cell_moment(vcat(SVector(z), Q), vl, vr), b + 1)
        end
    end
    return (MF0, MF1, MF2, MT0, MT1)
end

# One cyclotron harmonic: loop perp cells, precompute the parallel-moment
# t-polynomials once per cell, then a smooth p⊥ QuadGK whose integrand only
# evaluates those + Bessel weights.
function _grid_harmonic(n, fit::TensorSplineFit, ω, Ω, kz, a)
    kz0 = iszero(kz)
    c = kz0 ? 1 / (ω - n * Ω) : -1 / kz
    kq = fit.knots_perp
    acc = zero(SVector{6, ComplexF64})
    ζ = kz0 ? zero(ω) : (ω - n * Ω) / kz
    σ = sign(kz)
    for i in 1:(length(kq) - 1)
        wl, wr = kq[i], kq[i + 1]
        # function barrier: the two poly variants have different eltypes
        polys = kz0 ? _grid_parmoment_polys0(fit, i) : _grid_parmoment_polys(fit, i, ζ, σ)
        acc = acc .+ _grid_cell_integral(polys, wl, wr, c, n, a, ω, kz, n * Ω)
    end
    return acc
end

function _grid_cell_integral((MF0c, MF1c, MF2c, MT0c, MT1c), wl, wr, c, n, a, ω, kz, nΩ)
    integ = v -> begin
        t = v - wl
        M = (
            evalpoly(t, MF0c.data), evalpoly(t, MF1c.data), evalpoly(t, MF2c.data),
            evalpoly(t, MT0c.data), evalpoly(t, MT1c.data),
        )
        _In_block(M, c, _perp_Bessel_bilinear(n, a, v), v, ω, kz, nΩ)
    end
    # The Bessel weight J_n(a v) has v-wavelength ≈ π/a; adaptive QuadGK over a
    # cell spanning many wavelengths can't resolve it. Pre-split the cell at the
    # oscillation scale (≥2 panels/wavelength) and let QuadGK adapt within each.
    return first(_quadgk_osc(integ, wl, wr, a))
end

# Adaptive QuadGK with oscillation-scale breakpoints for the Bessel kernel.
@inline function _quadgk_osc(f, wl, wr, a)
    nb = ceil(Int, abs(a) * (wr - wl) / (π / 2))
    pts = nb <= 1 ? (wl, wr) : range(wl, wr; length = nb + 1)
    return QuadGK.quadgk(f, pts...; rtol = 1.0e-6, norm = NORM)
end
