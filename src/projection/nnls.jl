# Ref: MPDES `annlsSp`/`nnlsSp` residual-driven knot insertion
# On an 81×61 grid: ~0.02 s, same fit residual/χ accuracy as a full global Kronecker NNLS (~700× its cost).

"Separable two-pass nonneg B-spline fit (positivity)."
Base.@kwdef struct NonnegBSpline{O} <: GridFitMethod
    tol::Float64 = 1.0e-3
    maxknots_par::Int = typemax(Int)
    maxknots_perp::Int = typemax(Int)
end

order(::NonnegBSpline{O}) where {O} = O

# Clamped knot vector for breakpoints `t` (length m+1 -> m cells).
function _clamped_knots(t, p::Integer)
    [fill(t[1], p); t; fill(t[end], p)]
end

# Knots picked per axis (adaptive 1-D fit on the max-norm slice),
# then the control coeffs by two sets of small 1-D NNLS instead of one Kronecker solve.
function fit_grid(m::NonnegBSpline, vpar, vperp, F)
    @assert all(>=(0), F) "NonnegBSpline needs f₀ ≥ 0"
    deg = order(m)
    ref_col = argmax(j -> norm(view(F, :, j)), axes(F, 2))
    knots_par = select_knots_1d(vpar, F[:, ref_col]; tol=m.tol, order=deg,
        maxknots=min(m.maxknots_par, length(vpar)))
    ref_row = argmax(i -> norm(view(F, i, :)), axes(F, 1))
    knots_perp = select_knots_1d(vperp, F[ref_row, :]; tol=m.tol, order=deg,
        maxknots=min(m.maxknots_perp, length(vperp)))
    Bpar = _collocation_matrix(knots_par, deg, vpar)     # nvpar × nb_par
    Bperp = _collocation_matrix(knots_perp, deg, vperp)  # nvperp × nb_perp
    nb_par, nb_perp = size(Bpar, 2), size(Bperp, 2)
    # Pass 1 (v∥): perp columns as RHS → C1[p,j] = para control p at v⊥ = vperp[j].
    C1 = zeros(nb_par, length(vperp))
    for j in axes(F, 2)
        C1[:, j] .= nnls(Bpar, F[:, j])
    end
    # Pass 2 (v⊥): fit each para-control row along v⊥ → ctrl[p,:] (nonneg).
    ctrl = zeros(nb_par, nb_perp)
    for p in 1:nb_par
        ctrl[p, :] .= nnls(Bperp, C1[p, :])
    end
    return _ctrl_to_coeffs(knots_par, knots_perp, ctrl, deg)
end

# --- NNLS: min ||A c - b||_2  s.t. c >= 0 ---
# alg`:nnls`  beat `:fnnls`/`:pivot` for our small-n shapes.
nnls(A::AbstractMatrix, b::AbstractVector) = vec(nonneg_lsq(A, b; alg=:nnls))

# --- 1-D adaptive NNLS knot selection ---

# Evaluate B-spline with control coeffs `c` at `x` (used in B-spline→power conversion).
function _bspline_eval(breakpoints, p, c, x)
    knots = _clamped_knots(breakpoints, p)
    nb = length(knots) - p - 1
    N = _bspline_basis_all(knots, p, x, nb)
    sum(N[i] * c[i] for i in 1:nb)
end

# Convert nonneg-B-spline control coeffs to per-cell power-basis [a,b,c,d] via
# local Vandermonde solve at 4 equispaced nodes per cell (degree-3 exact: 4
# points determine a cubic exactly, no least-squares needed).
function _bspline_to_power(breakpoints, p, c)
    ncell = length(breakpoints) - 1
    coeffs = zeros(eltype(c), ncell, p + 1)
    for i in 1:ncell
        v0, v1 = breakpoints[i], breakpoints[i+1]
        h = v1 - v0
        s_nodes = range(0, h; length=p + 1)  # p+1=4 nodes, exact for a cubic
        fvals = [_bspline_eval(breakpoints, p, c, v0 + s) for s in s_nodes]
        V = [s^k for s in s_nodes, k in 0:p]
        coeffs[i, :] .= V \ fvals
    end
    coeffs
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
function select_knots_1d(v::AbstractVector, f::AbstractVector; tol::Real=1e-3, order::Integer=3, maxknots::Integer=length(v))
    @assert all(>=(0), f) "input f0 must be nonnegative"
    n = length(v)
    fnorm = norm(f)
    fnorm == 0 && return [v[1], v[end]]

    breakpoints = [v[1], v[end]]
    while true
        B = _collocation_matrix(breakpoints, order, v)
        resid = f .- B * nnls(B, f)
        ncells = length(breakpoints) - 1
        (norm(resid) / fnorm <= tol || ncells >= maxknots || ncells >= n - 1) && return breakpoints
        # bin squared residual into cells, split the worst one at its
        # data-index midpoint (mirrors annlsSp's newKnots cumulative-residual rule)
        res2 = resid .^ 2
        worst_cell, worst_val = 1, -Inf
        for i in 1:ncells
            lo, hi = breakpoints[i], breakpoints[i+1]
            s = sum(res2[(v .>= lo) .& (v .<= hi)])
            if s > worst_val
                worst_val = s
                worst_cell = i
            end
        end
        lo, hi = breakpoints[worst_cell], breakpoints[worst_cell+1]
        idxs = findall(x -> lo < x < hi, v)
        newpt = isempty(idxs) ? (lo + hi) / 2 : v[idxs[(length(idxs)+1)÷2]]
        insert!(breakpoints, worst_cell + 1, newpt)
    end
end

# Collocation matrix B[k,i] = B_{i,p}(x[k]) for nb = length(breakpoints)+p-1 basis functions.
function _collocation_matrix(breakpoints::AbstractVector, p::Integer, xs::AbstractVector)
    knots = _clamped_knots(breakpoints, p)
    nb = length(knots) - p - 1
    B = zeros(Float64, length(xs), nb)
    for (k, x) in enumerate(xs)
        B[k, :] .= _bspline_basis_all(knots, p, x, nb)
    end
    B
end

# Per-cell power-basis coefficients of a single tensor-product B-spline basis
# function B_i(v∥)*B_j(v⊥), as a (ncell_par, ncell_perp, 4, 4) array — built
# once per axis (outer product of the 1-D per-basis-function power tables) and
# reused for every control coefficient via linearity.
function _bspline_to_power_table(breakpoints, p)
    knots = _clamped_knots(breakpoints, p)
    nb = length(knots) - p - 1
    ncell = length(breakpoints) - 1
    # table[i][cell,:] = power coeffs of basis function i restricted to `cell`
    [_bspline_to_power(breakpoints, p, [k == i ? 1.0 : 0.0 for k in 1:nb]) for i in 1:nb]
end

# Tensor B-spline control coeffs `ctrl[p,q]` → per-cell power coeffs (degree `p`).
# Each basis B_p(v∥)B_q(v⊥) contributes its (precomputed) per-cell power table
# outer product, scaled by `ctrl[p,q]`. Shared by every NNLS fit method.
function _ctrl_to_coeffs(knots_par, knots_perp, ctrl, order)
    ncell_par, ncell_perp = length(knots_par) - 1, length(knots_perp) - 1
    tab_par = _bspline_to_power_table(knots_par, order)    # nb_par tables, each ncell_par x (order+1)
    tab_perp = _bspline_to_power_table(knots_perp, order)  # nb_perp tables, each ncell_perp x (order+1)
    nb_par, nb_perp = length(tab_par), length(tab_perp)
    N = order + 1

    acc = [zeros(N, N) for _ in 1:ncell_par, _ in 1:ncell_perp]
    for p in 1:nb_par, q in 1:nb_perp
        w = ctrl[p, q]
        w == 0 && continue
        Tp, Tq = tab_par[p], tab_perp[q]  # ncell_par x N, ncell_perp x N
        for icell in 1:ncell_par, jcell in 1:ncell_perp
            ap = @view Tp[icell, :]
            aq = @view Tq[jcell, :]
            any(!=(0), ap) || continue
            any(!=(0), aq) || continue
            acc[icell, jcell] .+= w .* (ap * aq')
        end
    end

    coeffs = [SMatrix{N,N,Float64}(acc[i, j]) for i in 1:ncell_par, j in 1:ncell_perp]
    TensorSplineFit(knots_par, knots_perp, coeffs, ctrl)
end



# Cox-de Boor recursion, single basis function B_{i,p} at x (0-indexed i internally via 1-based here).
function _bspline_basis_all(knots::AbstractVector{T}, p::Integer, x::Real, nb::Integer) where {T}
    m = length(knots) - 1
    N = zeros(promote_type(T, typeof(x)), m) # degree-0 stage buffer, length m (knots has m+1 entries)
    for i in 1:m
        N[i] = (knots[i] <= x < knots[i+1]) || (x == knots[end] && knots[i] < knots[i+1] && i == m - p) ? one(eltype(N)) : zero(eltype(N))
    end
    # handle x exactly at right end robustly: ensure partition of unity at x=knots[end]
    if x == knots[end]
        fill!(N, 0)
        N[m] = 1
    end
    for deg in 1:p
        Nnew = zeros(eltype(N), m)
        for i in 1:(m-deg)
            left = knots[i+deg] - knots[i]
            right = knots[i+deg+1] - knots[i+1]
            a = left > 0 ? (x - knots[i]) / left * N[i] : zero(eltype(N))
            b = right > 0 ? (knots[i+deg+1] - x) / right * N[i+1] : zero(eltype(N))
            Nnew[i] = a + b
        end
        N = Nnew
    end
    @view N[1:nb]
end
