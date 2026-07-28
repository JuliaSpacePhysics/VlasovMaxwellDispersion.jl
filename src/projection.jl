# `fit_grid` dispatches on the method; add a method by implementing one `fit_grid`
#
# Two methods, by trade-off (both yield TensorSplineFit):
#   NonnegBSpline (default) : SEPARABLE two-pass nonneg
#     fit. Pass 1 fits the v∥ B-spline with the v⊥ columns as multi-RHS; pass 2
#     fits the resulting control rows along v⊥. No Kronecker matrix.
#     + f₀≥0 (nonneg coeffs each pass), C², O(per-axis) small NNLS.
#   BicubicHermite — local C¹ Hermite interpolation, knots = grid, FD-estimated
#     ∂∥,∂⊥,∂∥∂⊥ at nodes; per-cell patch = M·Q·Mᵀ (no solve at all).
#     + O(ngrid), cheapest, exact at nodes, reproduces any tensor-cubic.

using LinearAlgebra: cholesky, Symmetric, norm

"""
    fit_grid(method, vperp, vpara, F) -> projection

Project gridded `F[i,j] = f₀(vperp[i], vpara[j])` onto a callable representation.
It may specialize `_grad2` and `moment` for analytic derivatives and integration. 
Assume `vperp`/`vpara` ascending, `vperp[1] ≥ 0`.
"""
function fit_grid end

abstract type GridFitMethod end

"""`domain(x) -> (; perp=(lo, hi), para=(lo, hi))`"""
function domain end

"""
    moment(projection, (mperp, mpara))

Compute `∫∫ vperp^mperp * vpara^mpara * projection(vperp, vpara)
dvpara dvperp` over `domain(projection)`.
"""
function moment(projection, powers::NTuple{2,Integer}; rtol=1.0e-9)
    mperp, mpara = powers
    min(mperp, mpara) >= 0 || throw(ArgumentError("moment powers must be nonnegative"))
    lims = domain(projection)
    return first(QuadGK.quadgk(lims.perp...; rtol) do v
        v^mperp * QuadGK.quadgk(u -> u^mpara * projection(v, u), lims.para...; rtol)[1]
    end)
end

include("projection/nnls.jl")
include("projection/bicubic_hermite.jl")

"""
    TensorSplineFit

Tensor spline fit on a `(v_perp, v_par)` grid. Per-`(i,j)` cell, `coeffs[i,j]` is
a `NP×NQ` matrix `c[A,B]` for
`f(v⊥,v∥) = Σ_{A,B} c[A,B] s⊥^{A-1} s∥^{B-1}`, `s⊥=v⊥-knots_perp[i]`, `s∥=v∥-knots_para[j]`.
Perpendicular degree `NP-1`, parallel degree `NQ-1`.
"""
struct TensorSplineFit{VP,VQ,T,NP,NQ,L}
    knots_perp::VP
    knots_para::VQ
    coeffs::Matrix{SMatrix{NP,NQ,T,L}}   # ncells_perp × ncells_par, each a per-cell power-coeff matrix
end

# the Landau contour pushes v∥ off-axis while the piecewise domain is partitioned on the real line.
@inline _cell(knots, x) = clamp(searchsortedlast(knots, real(x)), 1, length(knots) - 1)

domain(fit::TensorSplineFit) = (;
    perp=(first(fit.knots_perp), last(fit.knots_perp)),
    para=(first(fit.knots_para), last(fit.knots_para)),
)

@inline function _insupport(fit::TensorSplineFit, v, u)
    lims = domain(fit)
    return lims.perp[1] <= real(v) <= lims.perp[2] &&
           lims.para[1] <= real(u) <= lims.para[2]
end

# Spline value / derivatives at (v,u). v⊥ real; u∥ may be complex (parallel Landau
# continuation). coeffs[i,j][A,B] multiplies s⊥^{A-1} s∥^{B-1}.
# Out-of-support zero MUST match the in-support `_polyval2` element type (real for real
# (v,u), complex only when an arg is).
function (fit::TensorSplineFit)(v, u)
    _insupport(fit, v, u) || return zero(promote_type(typeof(float(v)), typeof(float(u)), eltype(eltype(fit.coeffs))))
    i, j = _cell(fit.knots_perp, v), _cell(fit.knots_para, u)
    return _polyval2(fit.coeffs[i, j], v - fit.knots_perp[i], u - fit.knots_para[j])
end

@inline function _polyval2(c, s, t)
    acc = zero(promote_type(typeof(s), typeof(t), eltype(c)))
    @inbounds for B in axes(c, 2), A in axes(c, 1)
        acc += c[A, B] * s^(A - 1) * t^(B - 1)
    end
    return acc
end

@inline function _grad2(fit::TensorSplineFit, v, u)
    _insupport(fit, v, u) ||
        (z=zero(promote_type(typeof(float(v)), typeof(float(u)), eltype(eltype(fit.coeffs)))); return (z, z))
    i, j = _cell(fit.knots_perp, v), _cell(fit.knots_para, u)
    return _dgradpolyval2(fit.coeffs[i, j], v - fit.knots_perp[i], u - fit.knots_para[j])
end

# Fused (∂s, ∂t) of the cell poly in one pass.
@inline function _dgradpolyval2(c, s, t)
    T = promote_type(typeof(s), typeof(t), eltype(c))
    ds = zero(T)
    dt = zero(T)
    tB1 = one(T)            # t^(B-1)
    tB2 = zero(T)           # t^(B-2); B=1 column has no ∂t
    @inbounds for B in axes(c, 2)
        sA1 = one(T)        # s^(A-1)
        sA2 = zero(T)       # s^(A-2); A=1 row has no ∂s
        for A in axes(c, 1)
            cAB = c[A, B]
            ds += (A - 1) * cAB * sA2 * tB1
            dt += (B - 1) * cAB * sA1 * tB2
            sA2, sA1 = sA1, sA1 * s
        end
        tB2, tB1 = tB1, tB1 * t
    end
    return (ds, dt)
end

@inline function _power_interval_moment(lo, hi, local_degree, power)
    h = hi - lo
    acc = zero(float(promote_type(typeof(lo), typeof(hi))))
    @inbounds for r in 0:power
        acc += binomial(power, r) * lo^(power - r) *
               h^(local_degree + r + 1) / (local_degree + r + 1)
    end
    return acc
end

function moment(fit::TensorSplineFit, powers::NTuple{2,Integer}; kw...)
    mperp, mpara = powers
    min(mperp, mpara) >= 0 || throw(ArgumentError("moment powers must be nonnegative"))
    kq, kp, coeffs = fit.knots_perp, fit.knots_para, fit.coeffs
    T = promote_type(eltype(eltype(coeffs)), eltype(kq), eltype(kp))
    acc = zero(T)
    @inbounds for i in axes(coeffs, 1), j in axes(coeffs, 2)
        cell = coeffs[i, j]
        for A in axes(cell, 1), B in axes(cell, 2)
            iq = _power_interval_moment(kq[i], kq[i+1], A - 1, mperp)
            ip = _power_interval_moment(kp[j], kp[j+1], B - 1, mpara)
            acc += cell[A, B] * iq * ip
        end
    end
    return acc
end
