# Exact per-k plan for CoupledVDF: the fixed-node analog of `LowRankPlan` with no surrogate.
# Cost is O(nperp·nmax) per ω vs LowRank's O(rank·nmax).
# Selected by `backend = FixedNodeEval()`.

struct CoupledPlan{D,T,S2,S5,S6}
    d::D                                # CoupledVDF (raw f₀, dgrad)
    Ω::T; kz::T; Pi2::T; invn::T
    ns::UnitRange{Int}
    U::T                                # p∥ half-width: far/near split
    lims::Tuple{T,T}                    # (p∥ lo, hi)
    un::Vector{T}; uw::Vector{T}        # fixed p∥ nodes/weights (near field)
    QP::Matrix{S2}                      # [l,j] (∂⊥f, ∂∥f) at (uₗ, vⱼ)
    nu::Matrix{S2}                      # [p,j] far scaled p∥ moments
    raw::Vector{S5}                     # [j] raw p∥ moments (q, uq, u²q, p, up)
    vn::Vector{T}; vw::Vector{T}        # fixed p⊥ nodes/weights
    b2::Matrix{S6}                      # [n,j] perp Bessel bilinears (ω-independent)
end

# marginal proxy for adaptive panel discovery: structure joint in (v,u) but invisible
# here is the horizon past which the fixed grid under-resolves and AdaptiveEval is needed.
_marg(g, lo, hi, np) = sum(x -> abs(g(x)), range(lo, hi, np))

function _build_coupled_plan(c::PreparedVDF{<:CoupledVDF}, s, k, b::FixedNodeEval)
    d = c.vdf
    Ω, kz = s.Omega, para(k)
    a = perp(k) / Ω
    nmax = nmax_bessel(a^2 * abs(c.cache.pperp2_mean) / 2)
    ns = (-nmax):nmax
    ulo, uhi = d.para
    T = typeof(float(uhi))
    U = max(abs(ulo), abs(uhi))
    xg, wg = QuadGK.gauss(b.order)
    np = b.nprobe

    un, uw = gl_nodes(_panels(u -> _marg(v -> d.f0(v, u), d.perp..., np), ulo, uhi, zero(U)), xg, wg, T)
    # perp panels refined to the Bessel half-wavelength π/(2|a|), like the LowRank plan
    vn, vw = gl_nodes(_panels(v -> _marg(u -> d.f0(v, u), ulo, uhi, np), d.perp..., abs(a)), xg, wg, T)
    L = length(un); Nv = length(vn)

    S2 = SVector{2,T}
    QP = [S2(d.dgrad(vn[j], un[l])) for l in 1:L, j in 1:Nv]
    # νₚ[j] = (1/U)∫g(vⱼ,u)(u/U)^{p-1}du, taken on the SAME nodes as the near field. Unlike
    # `LowRankPara`, no adaptive moment rule: the far branch cannot beat the node sum the near
    # branch already commits to (measured: identical χ), and νₚ only ever enters weighted by
    # (U/ζ)ᵖ ≤ 2⁻ᵖ. One matrix product replaces Nv adaptive quadratures.
    tw = Matrix{T}(undef, _LR_NMOM + 1, L)
    for l in 1:L
        t = un[l] / U; q = uw[l] / U
        for p in axes(tw, 1)
            tw[p, l] = q; q *= t
        end
    end
    nu = tw * QP
    # ∫uᵖ g du = U^{p+1}·νₚ ⇒ the moment recursion's raw moments come free
    raw = [SVector(U * nu[1, j][1], U^2 * nu[2, j][1], U^3 * nu[3, j][1],
        U * nu[1, j][2], U^2 * nu[2, j][2]) for j in 1:Nv]

    S6 = SVector{6,T}
    b2 = Matrix{S6}(undef, length(ns), Nv)
    tmp = Vector{S6}(undef, length(ns))
    for j in 1:Nv
        _perp_Bessel_bilinears!(tmp, a, vn[j])
        @views b2[:, j] .= tmp
    end
    return CoupledPlan(d, Ω, kz, s.Pi2, one(T) / c.cache.n, ns, U, (T(ulo), T(uhi)),
        un, uw, QP, nu, raw, vn, vw, b2)
end

# Landau–Cauchy (A_j, B_j) = ∫(∂⊥f, ∂∥f)(vⱼ,u)/(u−ζ)du for every perp node j at one ζ,
# σ=sign(k∥). Same peel/Neumann split as `_lr_cauchy!`, but the columns are perp nodes
# (each an exact f₀ slice) instead of surrogate rank factors; the crossed residue uses
# d.dgrad(vⱼ,ζ) directly. The node kernel wₗ=uwₗ/(uₗ−ζ) is shared across all columns.
function _coupled_cauchy_cols!(AB, w, cs, pl::CoupledPlan, ζ, σ)
    lo, hi = pl.lims; U = pl.U; d = pl.d
    Z = eltype(AB)
    if abs(ζ) > _LR_THETA * U
        # Neumann weights (U/ζ)ᵖ up front: Horner per column would serialize 41 complex
        # multiplies on the latency-bound critical path, and this loop is the far-field cost.
        invξ = U / ζ; q = one(ζ)
        @inbounds for i in eachindex(cs)
            cs[i] = q *= invξ
        end
        crossed = σ * imag(ζ) < 0 && lo < real(ζ) < hi
        @inbounds for j in eachindex(AB)
            m = zero(Z)
            @simd for i in eachindex(cs)
                m += cs[i] * pl.nu[i, j]
            end
            AB[j] = crossed ? (σ * _2πim) * SVector(d.dgrad(pl.vn[j], ζ)) - m : -m
        end
        return
    end
    W, w0, hit = _node_kernel!(w, pl.un, pl.uw, ζ, _NODE_BAND * U)
    x = real(ζ); h = 1.0e-4 * (1 + abs(x))
    lb = searchsortedfirst(pl.un, clamp(x, lo, hi))
    l1, l2 = clamp(lb - 1, 1, length(pl.un)), clamp(lb, 1, length(pl.un))
    # both pole terms depend on ζ alone — only the peel GATE varies per column, so the
    # complex log must not sit inside the loop (it was the near-field cost)
    lpy, lpn = _lpole_term(ζ, lo, hi, σ, true), _lpole_term(ζ, lo, hi, σ, false)
    @inbounds for j in eachindex(AB)
        v = pl.vn[j]
        φζ = SVector(d.dgrad(v, ζ))
        # removable value standing in for the nodes `_node_kernel!` dropped: ∂ᵤ of the numerator
        dφζ = hit ? (SVector(d.dgrad(v, x + h)) .- SVector(d.dgrad(v, x - h))) ./ 2h : φζ
        S = zero(Z)
        @simd for l in eachindex(w)
            S += w[l] * pl.QP[l, j]
        end
        # the peel gate only needs |φ| on the real axis beside ζ — read it off the node table
        scale = max.(abs.(pl.QP[l1, j]), abs.(pl.QP[l2, j]))
        AB[j] = _coupled_near(S, W, w0, φζ, dφζ, scale, lpy, lpn)
    end
    return
end

# Finish both near-field transforms from their shared node sum: subtract φ(ζ)·Σₗwₗ where the
# pole is peeled, restore the removable node value, add the analytic pole term.
# Kept out of line so `W`/`w0` arrive as arguments — captured loop variables would box.
@inline _coupled_near(S, W, w0, φζ, dφζ, scale, lpy, lpn) =
    map(S, φζ, dφζ, scale) do s, f, df, sc
        _peel(f, sc) ? s + w0 * df - f * W + f * lpy : s + f * lpn
    end

# uʲ/(u−ζ) = ζʲ/(u−ζ) + poly ⇒ the u¹,u² Cauchy moments follow from A0,B0 and the raw
# moments; assemble each harmonic with the shared coupled kernels (`_In_forms`/`_In_assemble`).
function (pl::CoupledPlan)(ω)
    ωc = complex(float(ω))
    kz = pl.kz; z0 = iszero(kz)
    σ = z0 ? one(kz) : sign(kz)
    ik = z0 ? zero(kz) : -1 / kz
    T = typeof(ωc / oneunit(kz))
    Nv = length(pl.vn)
    AB = Vector{SVector{2,T}}(undef, Nv)
    w = Vector{T}(undef, length(pl.un))
    cs = Vector{T}(undef, _LR_NMOM + 1)
    acc = zero(SVector{6,ComplexF64})
    @inbounds for (i, nh) in enumerate(pl.ns)
        nΩ = nh * pl.Ω
        ζ = z0 ? ωc : (ωc - nΩ) / kz
        z0 || _coupled_cauchy_cols!(AB, w, cs, pl, ζ, σ)
        for j in 1:Nv
            r = pl.raw[j]
            Δm = if z0
                r .* inv(ωc - nΩ)
            else
                A0, B0 = AB[j]
                A1 = r[1] + ζ * A0
                A2 = r[2] + ζ * A1
                B1 = r[4] + ζ * B0
                ik .* SVector(A0, A1, A2, B0, B1)
            end
            F = _In_forms(Δm, pl.vn[j], ωc, kz)
            acc += (2π * pl.vw[j]) .* _In_assemble(F, pl.b2[i, j], nΩ, ωc)
        end
    end
    return (pl.Pi2 * pl.invn / ωc^2) * _antisymmat(acc)
end

# FixedNodeEval builds the fixed-node plan (non-relativistic only); AdaptiveEval and the
# relativistic regime fall back to the per-ω adaptive `GenericKPlan`.
function plan_contribution(c::PreparedVDF{<:CoupledVDF}, s, k;
    backend::ChiBackend = FixedNodeEval(), closure = HarmonicSum(), kw...)
    (backend isa FixedNodeEval && regime(c) isa NonRelativistic) || return GenericKPlan(s, k, closure)
    return _build_coupled_plan(c, s, k, backend)
end
plan_contribution(d::CoupledVDF, s, k; closure = HarmonicSum(), kw...) =
    plan_contribution(prepare(d, closure), s, k; closure, kw...)
