# --- Parallel Hilbert primitive H∥ for a piecewise-polynomial g ---
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

Synthetic division gives `p(v) = q(v)(v−ζ) + p(ζ)`; `∫q` is a plain polynomial
integral and `p(ζ)·log((vr−ζ)/(vl−ζ))` carries the branch cut as one complex log
of the ratio — the branch-cut invariant (spec §1).
"""
@inline function cell_hilbert(coeffs, vl, vr, ζ)
    m = length(coeffs)
    T = complex(promote_type(eltype(coeffs), typeof(float(vl)), typeof(float(vr)), typeof(ζ)))
    # Synthetic division of p by (v−ζ): quotient q has degree m-2, ascending in
    # qhi..qlo built top-down; p(ζ) is the final remainder (Horner on ζ).
    pζ = zero(T) + coeffs[m]          # running Horner value = current quotient coeff
    poly = zero(T)                    # ∫ q over the cell, accumulated by Horner-in-v
    # Walk quotient coeffs from highest (q_{m-2}) down to q_0. After processing
    # quotient coeff of degree d (=q value before update), its monomial v^d
    # contributes (vr^{d+1}-vl^{d+1})/(d+1) to ∫q. We integrate as we descend.
    for k in (m - 1):-1:1
        d = k - 1                     # degree of this quotient coefficient = pζ
        poly += pζ * (vr^(d + 1) - vl^(d + 1)) / (d + 1)
        pζ = coeffs[k] + ζ * pζ       # Horner step → remainder/next quotient coeff
    end
    # pζ now holds p(ζ). One complex log of the ratio = branch-cut-safe.
    poly + pζ * log((vr - ζ) / (vl - ζ))
end

"""
    hilbert_pwpoly(coeffs, nodes, ζ) -> Complex

Parallel Hilbert integral `H∥[g](ζ) = ∫ g(v)/(v − ζ) dv` for a piecewise
polynomial `g`: `nodes` are the `N+1` cell boundaries (ascending), `coeffs[i]`
is the monomial-coefficient vector (ascending degree) of `g` on cell
`[nodes[i], nodes[i+1]]`. Sums `cell_hilbert` over cells.

Landau-causal and single-valued across `Im ζ → 0` by construction (each cell
uses the log-of-ratio form). For `ζ` inside the support the physical sheet is
`Im ζ > 0`; the result continues analytically to `Im ζ < 0`.
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
    cell_hilbert_landau(coeffs, vl, vr, ζ) -> Complex

MPDES-style Landau-causal continuation of `cell_hilbert`. For `Im ζ < 0` and
`Re ζ` inside the cell, add `2πi p(ζ)` to continue from the upper half-plane.
"""
function cell_hilbert_landau(coeffs, vl, vr, ζ)
    h = cell_hilbert(coeffs, vl, vr, ζ)
    imag(ζ) < 0 && _pole_in_cell(vl, vr, ζ) ? h + 2π * im * _polyval(coeffs, ζ) : h
end

"""
    hilbert_landau_pwpoly(coeffs, nodes, ζ) -> Complex

Piecewise-polynomial Cauchy integral on the Landau sheet. This is the MPDES
parallel-pole correction: direct cell integral plus `2πi` residue for lower-half
poles whose real part lies in a cell.
"""
function hilbert_landau_pwpoly(coeffs, nodes, ζ)
    s = cell_hilbert_landau(coeffs[1], nodes[1], nodes[2], ζ)
    @inbounds for i in 2:length(coeffs)
        s += cell_hilbert_landau(coeffs[i], nodes[i], nodes[i + 1], ζ)
    end
    s
end
