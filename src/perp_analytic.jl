# --- Perpendicular Bessel-moment primitive P⊥ for a piecewise-polynomial h ---
# P⊥ needs cell integrals  ∫ v⊥^{d+1} J_n1(a v⊥) J_n2(a v⊥) dv⊥  (the v⊥^1 is the
# cylindrical Jacobian, d = monomial degree of h on the cell, a = k⊥/Ω_s the
# uniform Bessel argument coefficient).

# Ref:  MPDES `intPer`/`intPerSeries` (~1648,1815 in external/MPDES/MPDES.m).
#
# Closed form: the Bessel-product power series (Schläfli's ₂F₃ form; Watson, *A
# Treatise on the Theory of Bessel Functions* §5.41)
#   J_n1(z)J_n2(z) = Σ_k (-1)^k (z/2)^{n1+n2+2k} Γ(n1+n2+2k+1)
#                        / (k! Γ(n1+k+1)Γ(n2+k+1)Γ(n1+n2+k+1))
# is termwise a power of z, so ∫ v⊥^{d+1}·(term_k) dv⊥ is elementary
#
# DEAD ENDS:
# MPDES's `intPerAsymptotic` is wrong by 10-100% even at z~20 for n1≠n2
# (a transcription bug in the vendored MATLAB — its n==1 branch reuses `x_r` for `x_l`)
# and no closed 2-term Lommel form exists for n1≠n2

# Float64 series degrades smoothly: ~1e-10 rel. error by a·max(|vl|,|vr|)~10,
# total garbage by ~20 (calibrated against brute quadrature, see self-test).
const _PERP_SERIES_F64_MAX_Z = 8.0

# The ₂F₃ series is exact and cheap for small z but ALTERNATES with terms peaking
# at ~(z/2)^z, so Float64 loses all precision to cancellation past z~8 — an
# artifact of expanding from v=0, not of the integral (over one cell the integrand
# barely oscillates, a·h ≲ 1).

"""
    besselprod_moment(vl, vr, d, n1, n2, a; tol=1e-12) -> Real

Robust scalar perp cell integral `∫_{vl}^{vr} v^{d+1} J_{n1}(av) J_{n2}(av) dv`.
`a == 0` closes exactly (`J_0(0)=1`, higher orders vanish).
For small argument `a·max(|vl|,|vr|) ≤ _PERP_SERIES_F64_MAX_Z`
the Bessel-product power series is exact and cheapest.

Beyond that the alternating series loses Float64 precision to
cancellation; rather than re-run it in (slow) `BigFloat`, evaluate the smooth
integrand directly with stable Gauss–Legendre (`_besselprod_moment_gl`) — similar accuracy.
"""
function besselprod_moment(vl, vr, d, n1, n2, a; tol=1e-12)
    a == 0 && return (n1 == 0 && n2 == 0) * (vr^(d + 2) - vl^(d + 2)) / (d + 2)
    z = abs(a) * max(abs(vl), abs(vr))
    if z <= _PERP_SERIES_F64_MAX_Z
        return besselprod_moment_series(float(vl), float(vr), d, n1, n2, float(a); tol)
    end
    return _besselprod_moment_gl(float(vl), float(vr), d, n1, n2, float(a))
end


# Bessel-product power series (Schläfli ₂F₃), summed as a 4-term recurrence on the term coefficient
function besselprod_moment_series(vl, vr, d, n1, n2, a::T; tol=1e-14, kmax=4000) where {T}
    c = (a / 2)^(n1 + n2) / (gamma(T(n1 + 1)) * gamma(T(n2 + 1)))   # k=0 term coeff
    p = d + n1 + n2 + 2                                              # power of v in k=0 term
    term = c * (vr^p - vl^p) / p
    s = term
    k = 0
    while k < kmax
        # ratio of successive product-series coefficients c_{k+1}/c_k, exact
        ratio = -(a / 2)^2 * (n1 + n2 + 2k + 1) * (n1 + n2 + 2k + 2) /
                ((k + 1) * (n1 + k + 1) * (n2 + k + 1) * (n1 + n2 + k + 1))
        c *= ratio
        k += 1
        p += 2
        term = c * (vr^p - vl^p) / p
        s += term
        k > 3 && abs(term) <= abs(s) * tol && break
    end
    s
end


# Gauss–Legendre nodes for the large-z branch. 
# The integrand itself is smooth and — over one CELL
# barely oscillates (the product J_{n1}J_{n2} has phase 2av, so a cell of
# width h turns over only a·h/π times; for the perp grids here a·h ≲ 0.5)
# Panels are sized so a·(panel width) ≤ 1, then GL-16 (exact to degree 31) is ample.
const _BESSELPROD_GL = QuadGK.gauss(16)

function _besselprod_moment_gl(vl, vr, d, n1, n2, a)
    x, w = _BESSELPROD_GL
    m = max(1, ceil(Int, abs(a) * (vr - vl)))      # panels: a·(width/m) ≤ 1
    hh = (vr - vl) / m
    s = 0.0
    @inbounds for p in 0:(m-1)
        c = vl + (p + 0.5) * hh
        for i in eachindex(x)
            v = c + (hh / 2) * x[i]
            s += w[i] * v^(d + 1) * besselj(n1, a * v) * besselj(n2, a * v)
        end
    end
    (hh / 2) * s
end


"""
    cell_bessel_moment(coeffs, vl, vr, n1, n2, a) -> Real

`∫_{vl}^{vr} v·p(v)·J_{n1}(av) J_{n2}(av) dv` for one cell, `p(v) = Σ_k
coeffs[k] v^{k-1}`. Sums `besselprod_moment` over the monomial terms.
"""
function cell_bessel_moment(coeffs, vl, vr, n1, n2, a)
    s = zero(float(promote_type(eltype(coeffs), typeof(a))))
    @inbounds for k in eachindex(coeffs)
        d = k - 1
        iszero(coeffs[k]) && continue
        s += coeffs[k] * besselprod_moment(vl, vr, d, n1, n2, a)
    end
    s
end

"""
    perp_pwpoly(coeffs, nodes, n1, n2, a) -> Real

Perpendicular Bessel moment `P⊥[h](n1,n2,a) = ∫ v h(v) J_{n1}(av) J_{n2}(av) dv`
for a piecewise polynomial `h`: `nodes` are the cell boundaries, `coeffs[i]` the monomial coefficients of `h` on cell `[nodes[i], nodes[i+1]]`.
"""
function perp_pwpoly(coeffs, nodes, n1, n2, a)
    s = cell_bessel_moment(coeffs[1], nodes[1], nodes[2], n1, n2, a)
    @inbounds for i in 2:length(coeffs)
        s += cell_bessel_moment(coeffs[i], nodes[i], nodes[i+1], n1, n2, a)
    end
    s
end
