# Grid -> piecewise-poly basis (spec §1a)
using LinearAlgebra: cholesky, Symmetric, norm

"""
    fit_grid(method, vpar, vperp, F) -> TensorSplineFit

Build a tensor spline `f₀ ≈ Σ coeffs[i,j][A,B] s∥^{A-1} s⊥^{B-1}` from gridded
`F[i,j]=f₀(vpar[i],vperp[j])``. `vpar`/`vperp` ascending, `vperp[1] ≥ 0`.
"""
function fit_grid end

# --- Pluggable grid→bicubic fit methods --------------------------------------
# `fit_grid` dispatches on the method; add a method by implementing one `fit_grid`.
#
# Two methods, by trade-off (both yield TensorSplineFit):
#   NonnegBSpline (default) : SEPARABLE two-pass nonneg
#     fit. Pass 1 fits the v∥ B-spline with the v⊥ columns as multi-RHS; pass 2
#     fits the resulting control rows along v⊥. No Kronecker matrix.
#     + f₀≥0 (nonneg coeffs each pass), C², O(per-axis) small NNLS.
#   BicubicHermite — local C¹ Hermite interpolation, knots = grid, FD-estimated
#     ∂∥,∂⊥,∂∥∂⊥ at nodes; per-cell patch = M·Q·Mᵀ (no solve at all).
#     + O(ngrid), cheapest, exact at nodes, reproduces any tensor-cubic.
#     − only C¹, NO positivity guard (can undershoot <0 in steep tails), no
#       error-driven refinement; accuracy floored by grid resolution (FD derivs).
abstract type GridFitMethod end

include("projection/nnls.jl")
include("projection/bicubic_hermite.jl")

# --- 2-D tensor product (v_par, v_perp) projection ---
"""
    TensorSplineFit{T,NP,NQ,L}

Tensor spline fit on a `(v_par, v_perp)` grid. Per-`(i,j)` cell, `coeffs[i,j]` is
a static `NP×NQ` matrix `c[A,B]` for
`f(v∥,v⊥) = Σ_{A,B} c[A,B] s∥^{A-1} s⊥^{B-1}`, `s∥=v∥-knots_par[i]`, `s⊥=v⊥-knots_perp[j]`.
Parallel degree `NP-1`, perpendicular degree `NQ-1`.
Consumed cell-by-cell by `hilbert_pwpoly` (parallel polynomial in `s∥`) and
`perp_analytic` (transpose role for `s⊥`).
"""
struct TensorSplineFit{T,NP,NQ,L}
    knots_par::Vector{T}
    knots_perp::Vector{T}
    coeffs::Matrix{SMatrix{NP,NQ,T,L}}   # ncells_par × ncells_perp, each a per-cell power-coeff matrix
    ctrl::Matrix{T}                      # nb_par × nb_perp nonneg tensor B-spline control coefficients
end

# the Landau contour pushes v∥ off-axis while the piecewise domain is partitioned on the real line.
@inline _cell(knots, x) = clamp(searchsortedlast(knots, real(x)), 1, length(knots) - 1)

# Spline value / derivatives at (u,v). u may be complex (parallel Landau
# continuation); v real. coeffs[i,j][A,B] multiplies s∥^{A-1} s⊥^{B-1}.
# Out-of-support zero MUST match the in-support `_polyval2` element type (real for real
# (u,v), complex only when an arg is).
function (fit::TensorSplineFit)(u, v)
    _insupport(fit, u, v) || return zero(promote_type(typeof(float(u)), typeof(float(v)), eltype(eltype(fit.coeffs))))
    i, j = _cell(fit.knots_par, u), _cell(fit.knots_perp, v)
    _polyval2(fit.coeffs[i, j], u - fit.knots_par[i], v - fit.knots_perp[j])
end

# Σ c[A,B] s^{A-1} t^{B-1} and its s-/t-derivatives; bounds unrolled from the cell type.
@inline function _polyval2(c, s, t)
    acc = zero(promote_type(typeof(s), typeof(t), eltype(c)))
    @inbounds for B in axes(c, 2), A in axes(c, 1)
        acc += c[A, B] * s^(A - 1) * t^(B - 1)
    end
    acc
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
    (ds, dt)
end
