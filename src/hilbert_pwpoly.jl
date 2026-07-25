# Parallel Hilbert primitive H∥ for a piecewise-polynomial g.
# H∥[g](ζ) = ∫ g(v)/(v − ζ) dv with g piecewise-poly on a v-grid. Each cell
# [v_i, v_{i+1}] with poly p(v)=Σ c_k v^{k-1} integrates in closed form: divide
# p(v) by (v−ζ) → p(v)=q(v)(v−ζ)+p(ζ), so
#   ∫ p/(v−ζ) = ∫ q  +  p(ζ)·log((v_{i+1}−ζ)/(v_i−ζ)).
# The direct Cauchy integral uses the principal log branch. The Landau-causal
# continuation adds the same 2πi in-cell lower-half correction as MPDES.

"""
    cell_hilbert(coeffs, vl, vr, ζ) -> Complex

Closed-form `∫_{vl}^{vr} p(v)/(v − ζ) dv` for one cell, `p(v) = Σ_k coeffs[k] v^{k-1}`
(monomial basis, ascending degree). Exact for all complex `ζ` (including off-axis
and, by limit, on the real axis away from `[vl,vr]`).
"""
@inline function cell_hilbert(coeffs, vl, vr, ζ)
    m = length(coeffs)
    T = complex(promote_type(eltype(coeffs), typeof(float(vl)), typeof(float(vr)), typeof(ζ)))
    # One descending pass does both halves of p = q(v−ζ) + p(ζ): the Horner value is the
    # degree-d quotient coefficient before its update, so ∫q accumulates as we descend.
    pζ = zero(T) + coeffs[m]
    poly = zero(T)
    for k in (m - 1):-1:1
        d = k - 1
        poly += pζ * (vr^(d + 1) - vl^(d + 1)) / (d + 1)
        pζ = coeffs[k] + ζ * pζ
    end
    # log OF THE RATIO, not a difference of logs: that is what keeps the branch cut safe.
    poly + pζ * log((vr - ζ) / (vl - ζ))
end

"""
    hilbert_pwpoly(coeffs, nodes, ζ) -> Complex

Parallel Hilbert integral `H∥[g](ζ) = ∫ g(v)/(v − ζ) dv` for a piecewise
polynomial `g`: `nodes` are the `N+1` cell boundaries (ascending), `coeffs[i]`
is the monomial-coefficient vector (ascending degree) of `g` on cell
`[nodes[i], nodes[i+1]]`.

Single-valued across `Im ζ → 0` because every cell uses the log-of-ratio form. For `ζ`
inside the support the physical sheet is `Im ζ > 0`, continued analytically below.
"""
function hilbert_pwpoly(coeffs, nodes, ζ)
    s = cell_hilbert(coeffs[1], nodes[1], nodes[2], ζ)
    @inbounds for i in 2:length(coeffs)
        s += cell_hilbert(coeffs[i], nodes[i], nodes[i + 1], ζ)
    end
    s
end

@inline function _polyval(coeffs, x)
    acc = zero(complex(promote_type(eltype(coeffs), typeof(x))))
    @inbounds for k in length(coeffs):-1:1
        acc = coeffs[k] + x * acc
    end
    acc
end

@inline _pole_in_cell(vl, vr, ζ) = real(ζ) > vl && real(ζ) <= vr

"""
    cell_hilbert_landau(coeffs, vl, vr, ζ, σ=1) -> Complex

MPDES-style Landau-causal continuation of `cell_hilbert`. `σ = sign(k∥)` orients the
contour: for `σ·Im ζ < 0` (⟺ Im ω < 0) with `Re ζ` inside the cell, add `σ·2πi p(ζ)`
to continue from the causal half-plane.
"""
function cell_hilbert_landau(coeffs, vl, vr, ζ, σ = 1)
    h = cell_hilbert(coeffs, vl, vr, ζ)
    σ * imag(ζ) < 0 && _pole_in_cell(vl, vr, ζ) ? h + σ * 2π * im * _polyval(coeffs, ζ) : h
end

"""
    hilbert_landau_pwpoly(coeffs, nodes, ζ, σ=1) -> Complex

[`hilbert_pwpoly`](@ref) on the Landau sheet: [`cell_hilbert_landau`](@ref) summed over cells.
"""
function hilbert_landau_pwpoly(coeffs, nodes, ζ, σ = 1)
    s = cell_hilbert_landau(coeffs[1], nodes[1], nodes[2], ζ, σ)
    @inbounds for i in 2:length(coeffs)
        s += cell_hilbert_landau(coeffs[i], nodes[i], nodes[i + 1], ζ, σ)
    end
    s
end
