# Ref: MPDES `annlsSp`/`nnlsSp` residual-driven knot insertion

"Separable two-pass nonneg B-spline fit (positivity)."
Base.@kwdef struct NonnegBSpline{O} <: GridFitMethod
    rtol::Float64 = 1.0e-3
    maxknots_para::Int = typemax(Int)
    maxknots_perp::Int = typemax(Int)
end

order(::NonnegBSpline{O}) where {O} = O

# Clamped knot vector for breakpoints `t` (length m+1 -> m cells).
function _clamped_knots(t, p::Integer)
    return [fill(t[1], p); t; fill(t[end], p)]
end

# Knots picked per axis (adaptive 1-D fit on the max-norm slice),
# then the control coeffs by two sets of small 1-D NNLS instead of one Kronecker solve.
# F[i,j] = f₀(x[i], y[j])
function fit_grid(m::NonnegBSpline, x, y, F; rtol = m.rtol)
    @assert all(>=(0), F) "NonnegBSpline needs f₀ ≥ 0"
    deg = order(m)
    ref_col = argmax(j -> norm(view(F, :, j)), axes(F, 2))
    knotsx = select_knots_1d(
        x, F[:, ref_col]; rtol, order = deg,
        maxknots = min(m.maxknots_perp, length(x))
    )
    ref_row = argmax(i -> norm(view(F, i, :)), axes(F, 1))
    knotsy = select_knots_1d(
        y, F[ref_row, :]; rtol, order = deg,
        maxknots = min(m.maxknots_para, length(y))
    )
    Bx = _collocation_matrix(knotsx, deg, x)
    By = _collocation_matrix(knotsy, deg, y)
    nbx, nby = size(Bx, 2), size(By, 2)
    # Pass 1 (x)
    C1 = zeros(nbx, length(y))
    for j in axes(F, 2)
        C1[:, j] .= nnls(Bx, F[:, j])
    end
    # Pass 2 (y): fit each x-control row along y → ctrl[p,:] (nonneg)
    ctrl = zeros(nbx, nby)
    for p in 1:nbx
        ctrl[p, :] .= nnls(By, C1[p, :])
    end
    return _ctrl_to_coeffs(knotsx, knotsy, ctrl, deg)
end

# --- NNLS: min ||A c - b||_2  s.t. c >= 0 ---
# alg`:nnls`  beat `:fnnls`/`:pivot` for our small-n shapes.
nnls(A::AbstractMatrix, b::AbstractVector) = vec(nonneg_lsq(A, b; alg = :nnls))


# Evaluate B-spline with control coeffs `c` at `x` (used in B-spline→power conversion).
function _bspline_eval(breakpoints, p, c, x)
    knots = _clamped_knots(breakpoints, p)
    nb = length(knots) - p - 1
    N = _bspline_basis_all(knots, p, x, nb)
    return sum(N[i] * c[i] for i in 1:nb)
end

# Convert nonneg-B-spline control coeffs to per-cell power-basis via
# local Vandermonde solve at 4 equispaced nodes per cell
function _bspline_to_power(breakpoints, p, c)
    ncell = length(breakpoints) - 1
    coeffs = zeros(eltype(c), ncell, p + 1)
    for i in 1:ncell
        v0, v1 = breakpoints[i], breakpoints[i + 1]
        h = v1 - v0
        s_nodes = range(0, h; length = p + 1)  # p+1=4 nodes, exact for a cubic
        fvals = [_bspline_eval(breakpoints, p, c, v0 + s) for s in s_nodes]
        V = [s^k for s in s_nodes, k in 0:p]
        coeffs[i, :] .= V \ fvals
    end
    return coeffs
end

"""
    select_knots_1d(v, f; tol=1e-3, order=3, maxknots=length(v)) -> Vector

Residual-driven breakpoint selection for one axis. Starts from a 2-breakpoint
(1-cell) degree-`order` spline and refines: each iteration solves NNLS for the
control coefficients (positivity by construction), then bisects the cell with
largest cumulative squared residual, until the relative L2 residual
`||f-f_sp||/||f|| < tol` or `maxknots` cells are placed. Returns the breakpoints
only — the control coefficients are refit downstream by the tensor two-pass
NNLS, so no power-basis conversion happens here.
"""
function select_knots_1d(v::AbstractVector, f::AbstractVector; rtol = 1.0e-3, order::Integer = 3, maxknots::Integer = length(v))
    @assert all(>=(0), f) "input f0 must be nonnegative"
    n = length(v)
    fnorm = norm(f)
    fnorm == 0 && return [v[1], v[end]]

    breakpoints = [v[1], v[end]]
    while true
        B = _collocation_matrix(breakpoints, order, v)
        resid = f .- B * nnls(B, f)
        ncells = length(breakpoints) - 1
        (norm(resid) / fnorm <= rtol || ncells >= maxknots || ncells >= n - 1) && return breakpoints
        # bin squared residual into cells, split the worst one at its data-index midpoint.
        # Only cells with an interior data point are splittable
        res2 = resid .^ 2
        worst_cell, worst_val = 0, -Inf
        for i in 1:ncells
            lo, hi = breakpoints[i], breakpoints[i + 1]
            any(x -> lo < x < hi, v) || continue
            s = sum(res2[(v .>= lo) .& (v .<= hi)])
            if s > worst_val
                worst_val = s
                worst_cell = i
            end
        end
        worst_cell == 0 && return breakpoints   # every data point is a breakpoint
        lo, hi = breakpoints[worst_cell], breakpoints[worst_cell + 1]
        idxs = findall(x -> lo < x < hi, v)
        insert!(breakpoints, worst_cell + 1, v[idxs[(length(idxs) + 1) ÷ 2]])
    end
    return
end

# Collocation matrix B[k,i] = B_{i,p}(x[k]) for nb = length(breakpoints)+p-1 basis functions.
function _collocation_matrix(breakpoints::AbstractVector, p::Integer, xs::AbstractVector)
    knots = _clamped_knots(breakpoints, p)
    nb = length(knots) - p - 1
    B = zeros(Float64, length(xs), nb)
    for (k, x) in enumerate(xs)
        B[k, :] .= _bspline_basis_all(knots, p, x, nb)
    end
    return B
end

# Per-cell power-basis coefficients of a single tensor-product B-spline basis
# built once per axis (outer product of the 1-D per-basis-function power tables) and
# reused for every control coefficient via linearity.
function _bspline_to_power_table(breakpoints, p)
    knots = _clamped_knots(breakpoints, p)
    nb = length(knots) - p - 1
    ncell = length(breakpoints) - 1
    # table[i][cell,:] = power coeffs of basis function i restricted to `cell`
    return [_bspline_to_power(breakpoints, p, [k == i ? 1.0 : 0.0 for k in 1:nb]) for i in 1:nb]
end

# Tensor B-spline control coeffs `ctrl[p,q]` → per-cell power coeffs (degree `order`).
# Each basis B_p(axis1)B_q(axis2) contributes its (precomputed) per-cell power
# table outer product, scaled by `ctrl[p,q]`. Shared by every NNLS fit method.
function _ctrl_to_coeffs(knots1, knots2, ctrl, order)
    ncell1, ncell2 = length(knots1) - 1, length(knots2) - 1
    tab1 = _bspline_to_power_table(knots1, order)    # nb1 tables, each ncell1 x (order+1)
    tab2 = _bspline_to_power_table(knots2, order)    # nb2 tables, each ncell2 x (order+1)
    nb1, nb2 = length(tab1), length(tab2)
    N = order + 1

    acc = [zeros(N, N) for _ in 1:ncell1, _ in 1:ncell2]
    for p in 1:nb1, q in 1:nb2
        w = ctrl[p, q]
        w == 0 && continue
        Tp, Tq = tab1[p], tab2[q]  # ncell1 x N, ncell2 x N
        for icell in 1:ncell1, jcell in 1:ncell2
            ap = @view Tp[icell, :]
            aq = @view Tq[jcell, :]
            any(!=(0), ap) || continue
            any(!=(0), aq) || continue
            acc[icell, jcell] .+= w .* (ap * aq')
        end
    end

    coeffs = [SMatrix{N, N, Float64}(acc[i, j]) for i in 1:ncell1, j in 1:ncell2]
    return TensorSplineFit(knots1, knots2, coeffs, ctrl)
end


# Cox-de Boor recursion, single basis function B_{i,p} at x (0-indexed i internally via 1-based here).
function _bspline_basis_all(knots::AbstractVector{T}, p::Integer, x::Real, nb::Integer) where {T}
    m = length(knots) - 1
    N = zeros(promote_type(T, typeof(x)), m) # degree-0 stage buffer, length m (knots has m+1 entries)
    for i in 1:m
        N[i] = (knots[i] <= x < knots[i + 1]) || (x == knots[end] && knots[i] < knots[i + 1] && i == m - p) ? one(eltype(N)) : zero(eltype(N))
    end
    # handle x exactly at right end robustly: ensure partition of unity at x=knots[end]
    if x == knots[end]
        fill!(N, 0)
        N[m] = 1
    end
    for deg in 1:p
        Nnew = zeros(eltype(N), m)
        for i in 1:(m - deg)
            left = knots[i + deg] - knots[i]
            right = knots[i + deg + 1] - knots[i + 1]
            a = left > 0 ? (x - knots[i]) / left * N[i] : zero(eltype(N))
            b = right > 0 ? (knots[i + deg + 1] - x) / right * N[i + 1] : zero(eltype(N))
            Nnew[i] = a + b
        end
        N = Nnew
    end
    return @view N[1:nb]
end
