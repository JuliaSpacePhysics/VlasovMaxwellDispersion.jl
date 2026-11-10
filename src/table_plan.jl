# Fixed-node per-k susceptibility plan for every non-relativistic f₀ without a closed-form parallel
# factor. f₀'s gradient is split into COLUMNS c, each a perp tensor weight per harmonic
# (P∂, PF)[n,c] times a parallel pair (g_c, h_c): MF moments come from g_c, MT from h_c (`_chi_mblock`).
#   SeparableVDF f⊥⊗f∥ : one column (f∥, f∥′)
#   LowRankVDF         : one per rank (bₛ, bₛ′)
#   CoupledVDF         : one per coarse perp node vₘ, (∂⊥f, ∂∥f)(vₘ,·); the Bessel oscillation is
#                        folded into P by product integration against the vₘ interpolant, so the
#                        column count follows f₀'s structure, not k⊥.
# Columns are summed at plan time into Qₙ(uₗ) = Σ_c P[n,c]·g_c(uₗ), which is ω-independent: an ω costs
# O(nharm·L) node sums — a pole near the nodes is peeled with the table's own p∥ interpolant — plus
# exact values f₀(·, ζ) only for crossed poles (Landau residue). Those come from separate RESIDUE
# columns: an interpolant validated on the real axis degrades at complex ζ (a coupled f₀'s v-structure
# sharpens there), so a coupled f₀ evaluates them exactly at the Bessel-resolving perp nodes.

const S6{T} = SVector{6, T}

# Panel discovery only needs to LOCATE structure; the fixed GL rule per panel supplies the accuracy.
const _PANEL_RTOL = 1.0e-8

# marginal proxy ∫|g| for panel discovery
_marg(g, lo, hi, np) = sum(x -> abs(g(x)), range(lo, hi, np))

# Adaptive panel edges for `f` on [lo,hi]; `a>0` additionally caps panels at the Bessel
# half-wavelength π/(2a), which adaptive quadrature on a non-oscillatory proxy cannot see.
function _panels(f, lo, hi, a, rtol = _PANEL_RTOL)
    segs = QuadGK.quadgk_segbuf(f, lo, hi; rtol)[3]
    e = sort!(unique!(vcat([s.a for s in segs], [s.b for s in segs])))
    return _refine(e, a)
end

function _refine(e, a)
    iszero(a) && return e
    out = [e[1]]
    for i in 2:length(e)
        m = max(1, ceil(Int, a * (e[i] - e[i - 1]) / (π / 2)))
        append!(out, range(e[i - 1], e[i]; length = m + 1)[2:end])
    end
    return out
end

# Interpolation, not quadrature, sizes the tables: a GL panel that integrates f₀ to rtol
# interpolates it only to ~√rtol. In p⊥ (coarse columns of a coupled f₀) that error reaches χ
# undamped once the Bessel weights oscillate within the panel; in p∥ it is the pole-peel accuracy.

# Barycentric weights of the nodes `xs`; an affine map of `xs` only rescales them all.
_baryweights(xs) = [inv(prod(xs[j] - xs[k] for k in eachindex(xs) if k != j)) for j in eachindex(xs)]

# Lagrange cardinals of the nodes `xs` at `x`.
function _cardinals(xs, x)
    i = findfirst(==(x), xs)
    isnothing(i) || return [j == i ? one(x) : zero(x) for j in eachindex(xs)]
    t = _baryweights(xs) ./ (x .- xs)
    return t ./ sum(t)
end

#
# Bisect panels until interpolating the vector-valued `F` from the GL nodes `xg` reproduces it
# between nodes to `tol·max|F|`.
function _interp_panels(F, pan, xg, tol)
    # between nodes, and at the panel ends: the node polynomial peaks there
    xt = [-1; [(xg[i] + xg[i + 1]) / 2 for i in 1:(length(xg) - 1)]; 1]
    scale = maximum(v -> maximum(abs, F(v)), range(first(pan), last(pan), 64))
    minhalf = 1.0e-6 * (last(pan) - first(pan))
    out = [first(pan)]
    stack = [(pan[i], pan[i + 1]) for i in (length(pan) - 1):-1:1]
    while !isempty(stack)
        a, b = pop!(stack)
        mid, half = (a + b) / 2, (b - a) / 2
        xs = mid .+ half .* xg
        Fx = F.(xs)
        err = maximum(xt) do t
            x = mid + half * t
            maximum(abs, sum(_cardinals(xs, x) .* Fx) - F(x))
        end
        if err > tol * scale && half > minhalf
            push!(stack, (mid, b), (a, mid))
        else
            push!(out, b)
        end
    end
    return out
end

# A pole beyond `_FAR·U` (U = p∥ half-width) takes the node sum's own Laurent series
# Σₗ uwₗQₗ/(uₗ−ζ) = −Σₚ Mₚ·(U/ζ)ᵖ/ζ, truncated at `_FAR^-_NMOM` ≈ 6e-14. The moments are of the
# NODES, not of f₀, so the series and the node sum are one function and the switch is seamless.
const _FAR = 2.0
const _NMOM = 44

"""
k-independent parallel side of a [`TablePlan`](@ref): fixed p∥ nodes, the column pairs
`(g_c, h_c)` on them, their moments, and `φ(c, ζ) -> (g_c(ζ), h_c(ζ))` continued off the real axis.
`upan` only has to locate the columns' structure; it is bisected until they interpolate to `tol`.
"""
struct ParaTable{T, F}
    φ::F
    lims::Tuple{T, T}
    upan::Vector{T}
    un::Vector{T}; uw::Vector{T}
    λ::Vector{T}                 # barycentric weights of one panel's nodes
    G::Matrix{SVector{2, T}}     # [l,c]
    mom::Matrix{SVector{2, T}}   # [p,c] Σₗ uwₗ·(uₗ/U)ᵖ⁻¹·(g_c, h_c)(uₗ)
end

function ParaTable(φ, ncol, upan, order, tol)
    T = float(eltype(upan))
    xg, wg = QuadGK.gauss(order)
    upan = _interp_panels(u -> [x for c in 1:ncol for x in φ(c, u)], upan, xg, tol)
    un, uw = gl_nodes(upan, xg, wg, T)
    G = [SVector{2, T}(φ(c, un[l])) for l in eachindex(un), c in 1:ncol]
    U = max(abs(first(upan)), abs(last(upan)))
    mom = [sum(l -> (uw[l] * (un[l] / U)^(p - 1)) * G[l, c], eachindex(un)) for p in 1:(_NMOM + 1), c in 1:ncol]
    return ParaTable(φ, (T(first(upan)), T(last(upan))), collect(T, upan), un, uw, _baryweights(T.(xg)), G, mom)
end

# Bessel bilinears Kₙ(v) for every n in ns at each perp node: K[i,j] = Kₙᵢ(vⱼ) (`_symmat` order).
function _bessel_tensors(a, ns, vs)
    T = float(promote_type(typeof(a), eltype(vs)))
    nmax = maximum(abs, ns)
    K = Matrix{S6{T}}(undef, length(ns), length(vs))
    Jv = Vector{T}(undef, nmax + 2)
    for (j, v) in enumerate(vs)
        besselj_ladder!(Jv, nmax + 1, a * v)
        for (i, n) in enumerate(ns)
            Jm, Jn, Jp = _jladder(Jv, n - 1), _jladder(Jv, n), _jladder(Jv, n + 1)
            b1, b2 = v * (Jm + Jp) / 2, v * (Jm - Jp) / 2
            K[i, j] = S6(b1^2, b1 * b2, b1 * Jn, b2^2, b2 * Jn, Jn^2)
        end
    end
    return K
end

# P∂[n,c] = Σⱼ 2πwⱼ αⱼ꜀ Kₙ(vⱼ),  PF[n,c] = Σⱼ 2πwⱼ vⱼ βⱼ꜀ Kₙ(vⱼ): perp quadrature of the column
# weights α (pairs with g_c, the ∂⊥ slice) and β (pairs with h_c, the ∂∥ slice).
function _perp_tensors(K, vs, ws, α, β)
    T = eltype(eltype(K))
    Km = Matrix(reinterpret(T, K))                   # 6nharm × nnodes, BLAS-ready
    P∂ = Km * ((2π .* ws) .* α)
    PF = Km * ((2π .* ws .* vs) .* β)
    return copy(reinterpret(S6{T}, P∂)), copy(reinterpret(S6{T}, PF))
end

"""
    TablePlan

Fixed-node per-`k` plan `ω ↦ χ` over a [`ParaTable`](@ref) with perp tensors `(P∂, PF)[n,c]`.
`residue = (φ, P∂, PF)` gives the columns contracted at near/crossed poles, `φ(j, ζ) -> (g_j, h_j)(ζ)`;
by default the table's own.
"""
struct TablePlan{PT, F, T}
    para::PT
    φ::F
    P∂::Matrix{S6{T}}; PF::Matrix{S6{T}}   # [n,j]: residue columns
    Q∂::Matrix{S6{T}}; QF::Matrix{S6{T}}   # [l,n]: Σ_c P[n,c]·(g_c, h_c)(uₗ)
    M∂::Matrix{S6{T}}; MF::Matrix{S6{T}}   # [p,n]: Σ_c P[n,c]·mom[p,c]
    ns::UnitRange{Int}
    Ω::T; kz::T; pref::T                   # pref = Π̃²/n
end

function TablePlan(t::ParaTable, P∂, PF, ns, Ω, kz, pref; residue = (t.φ, P∂, PF))
    Gg = first.(t.G); Gh = last.(t.G)
    Q∂ = permutedims(_contract(P∂, Gg)); QF = permutedims(_contract(PF, Gh))
    M∂ = permutedims(_contract(P∂, first.(t.mom))); MF = permutedims(_contract(PF, last.(t.mom)))
    return TablePlan(t, residue..., Q∂, QF, M∂, MF, ns, promote(Ω, kz, pref)...)
end

# out[n,l] = Σ_c P[n,c]·X[l,c]
function _contract(P::Matrix{S6{T}}, X) where {T}
    Pm = reshape(reinterpret(T, P), 6 * size(P, 1), size(P, 2))
    return copy(reinterpret(S6{T}, Pm * transpose(X)))
end

# A pole within a panel width of the nodes is under-resolved by the node sum and is peeled.
@inline function _near(upan, p, ζ, lo, hi)
    x = real(ζ)
    band = max(upan[p + 1] - upan[p], upan[min(p + 2, end)] - upan[p + 1], upan[p + 1] - upan[max(p - 1, 1)])
    return abs(complex(x - clamp(x, lo, hi), imag(ζ))) < band
end

# Interpolant of the panel nodes `rng` of table column Q[:,i] at `z`.
@inline function _interp(Q, i, rng, un, λ, z)
    num = zero(eltype(Q)) * zero(z); den = zero(z)
    @inbounds for (j, l) in enumerate(rng)
        d = z - un[l]
        iszero(d) && return Q[l, i] * one(z)
        c = λ[j] / d
        num += c * Q[l, i]; den += c
    end
    return num / den
end

# Its derivative at the node `l`.
@inline function _dinterp(Q, i, rng, un, λ, l)
    d = zero(eltype(Q)); λl = λ[l - first(rng) + 1]
    @inbounds for (j, m) in enumerate(rng)
        m == l || (d += (λ[j] / λl) * (Q[m, i] - Q[l, i]) / (un[l] - un[m]))
    end
    return d
end

# Σₗ uwₗ·Q[l,i]/(uₗ−ζ) for a far pole, from the node moments.
@inline function _far_sums(pl, i, ζ, U)
    ξ = U / ζ
    S∂ = zero(eltype(pl.M∂)) * zero(ζ); SF = S∂
    # terms until |ξ|ᴺ < 1e-13
    @inbounds for q in min(_NMOM, ceil(Int, -60 / log(abs2(ξ)))):-1:0
        S∂ = S∂ * ξ + pl.M∂[q + 1, i]; SF = SF * ξ + pl.MF[q + 1, i]
    end
    return -S∂ / ζ, -SF / ζ
end

# The same by the node sum; a pole near the nodes is peeled (adds the log, not the crossed residue).
@inline function _node_sums(pl, w, i, ζ, σ, crossed, U)
    t = pl.para; lo, hi = t.lims; x = real(ζ)
    W, w0, hit = _node_kernel!(w, t.un, t.uw, ζ, _NODE_BAND * U)
    S∂ = zero(eltype(pl.Q∂)) * zero(ζ); SF = S∂
    @inbounds for l in eachindex(w)
        S∂ += w[l] * pl.Q∂[l, i]; SF += w[l] * pl.QF[l, i]
    end
    p = clamp(searchsortedlast(t.upan, x), 1, length(t.upan) - 1)
    _near(t.upan, p, ζ, lo, hi) || return S∂, SF
    # Peel c = (panel interpolant)(ζ): Σₗwₗ(Qₗ − c) + c·log. An identity for any c, and
    # this c is consistent with the Qₗ it is subtracted from, so a node grazing ζ
    # (wₗ ~ 1/δ) amplifies nothing — an exact f₀(ζ) would differ from the table by its
    # interpolation error, times 1/δ.
    order = length(t.λ)
    rng = ((p - 1) * order + 1):(p * order)
    lg = _lpole_term(ζ, lo, hi, σ, true) - W
    crossed && (lg -= σ * _2πim)
    S∂ += lg * _interp(pl.Q∂, i, rng, t.un, t.λ, ζ)
    SF += lg * _interp(pl.QF, i, rng, t.un, t.λ, ζ)
    if hit   # removable value of the node `_node_kernel!` dropped
        l = argmin(m -> abs(t.un[m] - x), rng)
        S∂ += w0 * _dinterp(pl.Q∂, i, rng, t.un, t.λ, l)
        SF += w0 * _dinterp(pl.QF, i, rng, t.un, t.λ, l)
    end
    return S∂, SF
end

function (pl::TablePlan)(ω)
    ωc = complex(float(ω)); kz = pl.kz; t = pl.para
    z0 = iszero(kz)
    σ = z0 ? one(kz) : sign(kz); ik = z0 ? zero(kz) : -1 / kz
    CT = typeof(ωc)
    lo, hi = t.lims
    U = max(abs(lo), abs(hi))
    w = Vector{CT}(undef, size(pl.Q∂, 1))
    Φ = Vector{SVector{2, CT}}(undef, size(pl.P∂, 2))
    acc = zero(S6{CT})
    @inbounds for (i, n) in enumerate(pl.ns)
        nΩ = n * pl.Ω
        # raw moments ∫(g, ug, u²g) and ∫(h, uh)
        R∂0, R∂1, R∂2 = pl.M∂[1, i], U * pl.M∂[2, i], U^2 * pl.M∂[3, i]
        RF0, RF1 = pl.MF[1, i], U * pl.MF[2, i]
        if z0
            iΔ = inv(ωc - nΩ)
            MF0, MF1, MF2 = R∂0 * iΔ, R∂1 * iΔ, R∂2 * iΔ
            MT0, MT1 = RF0 * iΔ, RF1 * iΔ
        else
            ζ = (ωc - nΩ) / kz
            x = real(ζ)
            crossed = σ * imag(ζ) < 0 && lo < x < hi
            S∂, SF = abs2(ζ) > (_FAR * U)^2 ? _far_sums(pl, i, ζ, U) : _node_sums(pl, w, i, ζ, σ, crossed, U)
            if crossed
                φ∂, φF = _column_values!(Φ, pl, i, ζ)
                S∂ += (σ * _2πim) * φ∂; SF += (σ * _2πim) * φF
            end
            # uʲ/(u−ζ) = ζʲ/(u−ζ) + poly ⇒ higher moments from the raw ones
            A1 = R∂0 + ζ * S∂
            MF0, MF1, MF2 = ik * S∂, ik * A1, ik * (R∂1 + ζ * A1)
            MT0, MT1 = ik * SF, ik * (RF0 + ζ * SF)
        end
        wF0 = ωc * MF0 - kz * MF1; wF1 = ωc * MF1 - kz * MF2
        acc += S6(
            wF0[1] + kz * MT0[1], im * (wF0[2] + kz * MT0[2]), wF1[3] + kz * MT1[3],
            wF0[4] + kz * MT0[4], im * (wF1[5] + kz * MT1[5]), nΩ * MF2[6] + (ωc - nΩ) * MT1[6],
        )
    end
    return (pl.pref / ωc^2) * _antisymmat(acc)
end

# Σⱼ P[n,j]·(g_j, h_j)(z) from exact residue-column values.
@inline function _column_values!(Φ, pl, i, z)
    @inbounds for c in eachindex(Φ)
        Φ[c] = SVector(pl.φ(c, z))
    end
    φ∂ = zero(eltype(pl.Q∂)) * zero(eltype(eltype(Φ))); φF = φ∂
    @inbounds for c in eachindex(Φ)
        φ∂ += pl.P∂[i, c] * Φ[c][1]; φF += pl.PF[i, c] * Φ[c][2]
    end
    return φ∂, φF
end
