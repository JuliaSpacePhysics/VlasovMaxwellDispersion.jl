# Exact fixed-node plan for a non-relativistic CoupledVDF (`backend = FixedNodeEval()`), a
# `TablePlan` whose columns are coarse perp nodes vₘ resolving f₀. The columns and their p∥ table are
# k-independent (built once in `prepare`); per k only the Bessel weights are integrated against the
# vₘ interpolant.

struct CoupledTable{T, P}
    backend::FixedNodeEval
    vpan::Vector{T}     # coarse perp panels, sized for interpolation
    vc::Vector{T}       # their GL nodes: one column each
    para::P             # ParaTable of (∂⊥f, ∂∥f)(vₘ, ·)
end

function CoupledTable(d::CoupledVDF, b::FixedNodeEval)
    ulo, uhi = d.para
    T = typeof(float(uhi))
    xg, wg = QuadGK.gauss(b.order)
    np = b.nprobe
    upan = _panels(u -> _marg(v -> d.f0(v, u), d.perp..., np), ulo, uhi, zero(T))
    un, _ = gl_nodes(upan, xg, wg, T)
    vpan = _panels(v -> _marg(u -> d.f0(v, u), ulo, uhi, np), d.perp..., zero(T))
    vpan = _interp_panels(v -> [x for u in un for x in d.dgrad(v, u)], vpan, xg, b.tol)
    vc, _ = gl_nodes(vpan, xg, wg, T)
    para = ParaTable((c, ζ) -> d.dgrad(vc[c], ζ), length(vc), upan, b.order, b.tol)
    return CoupledTable(b, vpan, vc, para)
end

# Bessel-resolving perp nodes (vⱼ, wⱼ) and interpolation weights α[j,m] = ℓₘ(vⱼ) carrying them onto
# the coarse columns; a panel already within the Bessel half-wavelength keeps its own nodes (α = I).
function _product_weights(t::CoupledTable, a)
    order = t.backend.order
    xg, wg = QuadGK.gauss(order)
    T = eltype(t.vc)
    vs, ws = T[], T[]
    blocks = Tuple{UnitRange{Int}, UnitRange{Int}, Matrix{T}}[]
    for p in 1:(length(t.vpan) - 1)
        cols = ((p - 1) * order + 1):(p * order)
        sub = _refine(t.vpan[p:(p + 1)], a)
        v, w = gl_nodes(sub, xg, wg, T)
        ℓ = length(sub) == 2 ? Matrix{T}(I, order, order) : stack(x -> _cardinals(t.vc[cols], x), v; dims = 1)
        push!(blocks, (length(vs) .+ eachindex(v), cols, ℓ))
        append!(vs, v); append!(ws, w)
    end
    α = zeros(T, length(vs), length(t.vc))
    for (rows, cols, ℓ) in blocks
        α[rows, cols] .= ℓ
    end
    return vs, ws, α
end

# FixedNodeEval builds the fixed-node plan (non-relativistic only); AdaptiveEval and the
# relativistic regime fall back to the per-ω adaptive `GenericKPlan`.
function plan_contribution(c::PreparedVDF{<:CoupledVDF}, s, k;
    backend::ChiBackend = FixedNodeEval(), kw...)
    (backend isa FixedNodeEval && regime(c) isa NonRelativistic) || return GenericKPlan(s, k, _closure(backend))
    t = get(c.cache, :table, nothing)
    (isnothing(t) || t.backend != backend) && (t = CoupledTable(c.vdf, backend))
    Ω = s.Omega
    a = k.k_perp / Ω
    nmax = nmax_bessel(a^2 * abs(c.cache.pperp2_mean) / 2)
    ns = (-nmax):nmax
    vs, ws, α = _product_weights(t, abs(a))
    K = _bessel_tensors(a, ns, vs)
    P∂, PF = _perp_tensors(K, vs, ws, α, α)
    d = c.vdf
    residue = ((j, ζ) -> d.dgrad(vs[j], ζ), (2π .* ws') .* K, (2π .* ws' .* vs') .* K)
    return TablePlan(t.para, P∂, PF, ns, Ω, k.k_para, s.Pi2 / c.cache.n; residue)
end
plan_contribution(d::CoupledVDF, s, k; backend = FixedNodeEval(), kw...) =
    plan_contribution(prepare(d; backend), s, k; backend, kw...)
