using RationalFunctionApproximation: Barycentric, poles as aaa_poles

"""
    AAA(; n=(20, 16), tol=1e-13, max_degree=150, stagnation=5)

Derivative-free rational-fit global solver: `aaa`-fit `1/det(ω̃²𝒟)` on an
`n = (nRe, nIm)` grid over the ω window; the fit's poles are the det's zeros.
Fitting the deflated det ([`wave_dispersion_tensor`](@ref)) rather than raw
`det 𝒟` matters: the raw det's `ω=0` pole is a formulation artifact that
screens nearby roots from the fit's far field — exactly the low-frequency
(Alfvén/EMIC) roots a wide ω box must keep.
"""
Base.@kwdef struct AAA
    n::Tuple{Int, Int} = (20, 16)
    tol::Float64 = 1.0e-13
    max_degree::Int = 150
    stagnation::Int = 5
end

# `keep` drops samples an approximate VDF cannot reliably resolve
# fitting them spends the fit's degree on the surrogate's own continuation error
function discover(alg::AAA, f, region; keep = Returns(true))
    ll, ur = region
    Z = [
        complex(x, y)
            for x in range(real(ll), real(ur); length = alg.n[1])
            for y in range(imag(ll), imag(ur); length = alg.n[2])
            if keep(complex(x, y))
    ]
    nev = length(Z)
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
    fit, converged = _aaa(F, Z; alg.tol, alg.max_degree, alg.stagnation)
    zs = filter!(z -> _in_box(region, z), aaa_poles(fit))
    return zs, nev, converged
end

# Greedy AAA (Nakatsukasa et al.): interpolate at the worst sample, take the barycentric
# weights from the null vector of the Loewner matrix over the remaining (test) samples.
#
# The weights come from the SVD of the Loewner matrix's R factor, NOT of the matrix itself:
# `svd` would build the m×n left factor that nothing here consumes, and QR + an n×n SVD gives
# the identical last right-singular vector (to 1e-15) for half the cost. That halves the
# survey's dominant term, since this fit runs once per wavevector.
function _aaa(y, z; tol, max_degree, stagnation)
    m = length(z)
    fmax = maximum(abs, y)
    T = promote_type(eltype(z), eltype(y))
    # A node's own Loewner row is dropped, so the fit can never outrank the surviving samples.
    N = min(max_degree + 1, (m + 1) >> 1)
    C = Matrix{T}(undef, m, N)      # Cauchy 1/(zᵢ − σₖ)
    L = Matrix{T}(undef, m, N)      # Loewner (yᵢ − fσₖ)·C
    A = Matrix{T}(undef, m, N)      # active-row copy, consumed by qr!
    num = Vector{T}(undef, m); den = Vector{T}(undef, m)
    rows = Vector{Int}(undef, m)
    test = trues(m)
    _, i0 = findmin(abs, y)         # first node: same seed as `approximate`
    σ = [i0]; test[i0] = false
    errs = Float64[]
    hist = Vector{Tuple{Vector{Int}, Vector{T}}}()
    for n in 1:N
        j = σ[n]
        @inbounds for i in 1:m
            δ = z[i] - z[j]
            C[i, n] = iszero(δ) ? 1 / eps() : 1 / δ   # δ=0 only on duplicate samples
            L[i, n] = (y[i] - y[j]) * C[i, n]
        end
        nt = 0
        @inbounds for i in 1:m
            test[i] && (nt += 1; rows[nt] = i)
        end
        nt == 0 && break
        idx = view(rows, 1:nt)
        Av = view(A, 1:nt, 1:n)
        copyto!(Av, view(L, idx, 1:n))
        w = svd!(qr!(Av).R).V[:, end]

        Cv = view(C, idx, 1:n)
        mul!(view(num, 1:nt), Cv, w .* view(y, σ))
        mul!(view(den, 1:nt), Cv, w)
        err, imax = -Inf, 1
        @inbounds for i in 1:nt
            e = abs(y[idx[i]] - num[i] / den[i])
            e > err && ((err, imax) = (e, i))
        end
        push!(errs, err)
        push!(hist, (copy(σ), copy(w)))

        s = _quitting_check(errs, stagnation, tol, fmax, N)
        if s != 0
            s > 0 && ((σ, w) = hist[s])                # unconverged: fall back on the best fit
            return Barycentric(z[σ], y[σ], w), s < 0
        end

        k = idx[imax]
        push!(σ, k); test[k] = false
    end
    # Sample budget spent without ever passing the quitting check: also unconverged.
    s = argmin(errs)
    σ, w = hist[s]
    return Barycentric(z[σ], y[σ], w), false
end

# `RationalFunctionApproximation.quitting_check`, reproduced so the fit stops exactly where
# `approximate` did: -1 converged (keep the last fit), 0 continue, n rewind to fit n.
function _quitting_check(errs, stagnation, tol, fmax, max_iter)
    n = length(errs)
    errs[end] <= tol * fmax && return -1
    stagnant = false
    min_err, min_k = findmin(errs)
    if n >= stagnation + 5 && min_err < fmax / 100
        plateau = 5min_err
        stagnant = all(plateau < e for e in last(errs, stagnation)) || (min_k < n - 2stagnation)
    end
    return (n >= max_iter || stagnant) ? argmin(errs) : 0
end

# Deflation pins a structural zero at ω=0 for kinetic species (ω̃²χ → 0). No
# value-based test rejects it (the raw residual also vanishes there), so gate
# geometrically: fit/polish park it at ~0, genuine low-ω roots sit far above.
_origin_gate(::AAA, diag) = 1.0e-6 * diag
