using RationalFunctionApproximation: approximate, poles as aaa_poles

"""
    AAA(; n=(20, 16), kw...)

Derivative-free rational-fit global solver: `aaa`-fit `1/det(ω̃²𝒟)` on an
`n = (nRe, nIm)` grid over the ω window; the fit's poles are the det's zeros.
Fitting the deflated det ([`wave_dispersion_tensor`](@ref)) rather than raw
`det 𝒟` matters: the raw det's `ω=0` pole is a formulation artifact that
screens nearby roots from the fit's far field — exactly the low-frequency
(Alfvén/EMIC) roots a wide ω box must keep. A fit saturating `max_degree`
flags the solution `:Partial`; shrink the window.
"""
Base.@kwdef struct AAA{K}
    n::Tuple{Int,Int} = (20, 16)
    kw::K = (; stagnation=10)
end

function discover(alg::AAA, f, region)
    ll, ur = region
    Z = [
        complex(x, y)
        for x in range(real(ll), real(ur); length=alg.n[1])
        for y in range(imag(ll), imag(ur); length=alg.n[2])
    ]
    F = map(f, Z)
    nvalid = 0
    @inbounds for v in F
        nvalid += isfinite(v) && !iszero(v)
    end
    if nvalid == length(F)
        map!(inv, F, F)
    else
        z = similar(Z, nvalid)
        y = similar(F, nvalid)
        i = 0
        @inbounds for j in eachindex(F, Z)
            v = F[j]
            (isfinite(v) && !iszero(v)) || continue
            i += 1
            z[i], y[i] = Z[j], inv(v)
        end
        Z, F = z, y
    end
    fit = approximate(F, Z; alg.kw...)
    zs = filter!(z -> _in_box(region, z), aaa_poles(fit))
    return zs, prod(alg.n)
end

# Deflation pins a structural zero at ω=0 for kinetic species (ω̃²χ → 0). No
# value-based test rejects it (the raw residual also vanishes there), so gate
# geometrically: fit/polish park it at ~0, genuine low-ω roots sit far above.
_origin_gate(::AAA, diag) = 1.0e-6 * diag
