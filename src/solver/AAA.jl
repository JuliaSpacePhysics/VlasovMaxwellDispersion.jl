using RationalFunctionApproximation: aaa, poles as aaa_poles, degree as aaa_degree

"""
    AAA(; n=(20, 16), tol=1e-13, max_degree=150)

Derivative-free rational-fit global solver: `aaa`-fit `1/det(ω̃²𝒟)` on an
`n = (nRe, nIm)` grid over the ω window; the fit's poles are the det's zeros.
Fitting the deflated det ([`wave_dispersion_tensor`](@ref)) rather than raw
`det 𝒟` matters: the raw det's `ω=0` pole is a formulation artifact that
screens nearby roots from the fit's far field — exactly the low-frequency
(Alfvén/EMIC) roots a wide ω box must keep. A fit saturating `max_degree`
flags the solution `:Partial`; shrink the window.
"""
Base.@kwdef struct AAA
    n::Tuple{Int, Int} = (20, 16)
    tol::Float64 = 1.0e-13
    max_degree::Int = 150
end

function _slice_zeros(alg::AAA, f, region)
    ll, ur = region
    Z = [
        complex(x, y)
            for x in range(real(ll), real(ur); length = alg.n[1])
            for y in range(imag(ll), imag(ur); length = alg.n[2])
    ]
    F = map(f, Z)
    ok = map(v -> isfinite(v) && !iszero(v), F)   # 1/f samples need finite f ≠ 0
    fit = aaa(inv.(F[ok]), Z[ok]; alg.max_degree, alg.tol)
    zs = filter(z -> _in_box(region, z), aaa_poles(fit))
    return zs, aaa_degree(fit) >= alg.max_degree
end

# Deflation pins a structural zero at ω=0 for kinetic species (ω̃²χ → 0). No
# value-based test rejects it (the raw residual also vanishes there), so gate
# geometrically: fit/polish park it at ~0, genuine low-ω roots sit far above.
_origin_gate(::AAA, diag) = 1.0e-6 * diag
