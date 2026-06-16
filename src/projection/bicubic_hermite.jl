"Local bicubic Hermite interpolation on the grid (O(N), C¹, no positivity guard)."
struct BicubicHermite <: GridFitMethod end

const _HERMITE_M = @SMatrix [1.0 0 0 0; 0 0 1 0; -3 3 -2 -1; 2 -2 1 1]

# Non-uniform symmetric difference (one-sided at the ends): dy/dx at node i.
@inline function _dcentral(x, y, i)
    i == 1 && return (y[2] - y[1]) / (x[2] - x[1])
    i == lastindex(x) && return (y[i] - y[i-1]) / (x[i] - x[i-1])
    (y[i+1] - y[i-1]) / (x[i+1] - x[i-1])
end

# Per-cell bicubic Hermite patch → power coeffs[A,B] for s∥^{A-1} s⊥^{B-1}.
# Q packs corner data (value / ∂⊥·hy / ∂∥·hx / ∂∥∂⊥·hx·hy); M·Q·Mᵀ is the
# normalized-coord coeff matrix, then de-normalized by hx^{A-1}hy^{B-1}.
@inline function _hermite_patch(f, fx, fy, fxy, hx, hy)
    Q = @SMatrix [f[1, 1] f[1, 2] fy[1, 1]*hy fy[1, 2]*hy
        f[2, 1] f[2, 2] fy[2, 1]*hy fy[2, 2]*hy
        fx[1, 1]*hx fx[1, 2]*hx fxy[1, 1]*hx*hy fxy[1, 2]*hx*hy
        fx[2, 1]*hx fx[2, 2]*hx fxy[2, 1]*hx*hy fxy[2, 2]*hx*hy]
    A = _HERMITE_M * Q * _HERMITE_M'
    SMatrix{4,4,Float64}(A[a, b] / (hx^(a - 1) * hy^(b - 1)) for a in 1:4, b in 1:4)
end

function fit_grid(::BicubicHermite, vpar, vperp, F)
    vp, vq = collect(float.(vpar)), collect(float.(vperp))
    np, nq = length(vp), length(vq)
    @assert size(F) == (np, nq)
    Fx = [_dcentral(vp, view(F, :, j), i) for i in 1:np, j in 1:nq]
    Fy = [_dcentral(vq, view(F, i, :), j) for i in 1:np, j in 1:nq]
    Fxy = [_dcentral(vp, view(Fy, :, j), i) for i in 1:np, j in 1:nq]  # ∂∥ of ∂⊥
    coeffs = Matrix{SMatrix{4,4,Float64,16}}(undef, np - 1, nq - 1)
    for i in 1:(np-1), j in 1:(nq-1)
        hx, hy = vp[i+1] - vp[i], vq[j+1] - vq[j]
        f = SMatrix{2,2}(F[i, j], F[i+1, j], F[i, j+1], F[i+1, j+1])
        fx = SMatrix{2,2}(Fx[i, j], Fx[i+1, j], Fx[i, j+1], Fx[i+1, j+1])
        fy = SMatrix{2,2}(Fy[i, j], Fy[i+1, j], Fy[i, j+1], Fy[i+1, j+1])
        fxy = SMatrix{2,2}(Fxy[i, j], Fxy[i+1, j], Fxy[i, j+1], Fxy[i+1, j+1])
        coeffs[i, j] = _hermite_patch(f, fx, fy, fxy, hx, hy)
    end
    # ctrl unused downstream for this method; keep the nodal grid as metadata.
    TensorSplineFit(vp, vq, coeffs, Matrix{Float64}(F))
end
