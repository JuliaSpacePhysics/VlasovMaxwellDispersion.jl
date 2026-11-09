# Motivation: At fixed k the ONLY ω-dependence of expensive integral in χ comes from the pole ζₙ=(ω−nΩ)/k∥;
# the perp Bessel is ω-independent. For a coupled f₀ the inner Landau integral still
# depends on p⊥, so perp quadrature cannot leave the ω loop.
# The separable surrogate makes the perp Bessel moments ω-independent tensors built once per k.

"""
    LowRankVDF(f0; para, perp, dgrad=nothing, rtol=1e-8, rmax=40, probe=(200, 401))
    LowRankVDF(d::CoupledVDF; kw...)

Separable surrogate `f₀ ≈ Σₛ ãₛ(p⊥)·bₛ(p∥)` of a coupled gyrotropic `f₀(p⊥,p∥)`, built by 
adaptive cross approximation to relative tolerance `rtol` (rank capped by `rmax`, 
pivots searched on a `probe = (nperp, npara)` grid).

Drop-in for [`CoupledVDF`](@ref): Landau residues need `rank` slices per crossed harmonic
instead of every Bessel-resolving perp node, so it gains most at large `k⊥`. Coupled VDFs are
typically very low rank — a bi-kappa and a spherical shell are both ≈ 10.

`f0` must be evaluable at complex `p∥`: each `bₛ(u)=f₀(vₛ,u)` is a literal `f₀` slice, so unlike a
fitted surrogate ([`SplineVDF`](@ref)), which exists only on the real axis, it can be continued.

$(RANGE_DOC)

# Accuracy

The cross is fitted on the REAL axis, so `rtol` bounds χ only where the p∥ integral stays there:
real ω, growing modes, `k∥=0`.

A CROSSED Landau pole is different — it picks up the residue `f₀(p⊥,ζ)` at complex
`ζ=(ω−nΩ)/k∥`, where the `ãₛ` are extrapolating. That error grows with `|Im ω / k∥|`, is NOT
bounded by `rtol`, and eventually makes the surrogate's det gain exact zeros that are not modes of
`f₀`. Surveys drop those via [`trusted`](@ref); [`trust_error`](@ref) measures the horizon.
"""
struct LowRankVDF{F, G, H, T, C} <: SolverVDF
    f0::F
    dgrad::G
    bd::H                  # (v,u) ↦ (f₀, ∂∥f₀) for the ω hot path
    vp::Vector{T}          # perp pivots  → parallel factors bₛ(u) = f₀(vp[s], u)
    up::Vector{T}          # para  pivots → perp factors  ãₛ(v) = Σᵣ f₀(v,up[r])·M[r,s]
    M::Matrix{T}
    para::Tuple{T, T}
    perp::Tuple{T, T}
    cache::C               # all k-independent work, done once at construction
end

regime(::LowRankVDF) = NonRelativistic()

"""Separation rank of the surrogate: the number of `ãₛ⊗bₛ` terms the cross kept."""
LinearAlgebra.rank(d::LowRankVDF) = length(d.vp)

function LowRankVDF(f0; para, perp, dgrad = nothing, rtol = 1.0e-8, rmax = 40, probe = (200, 401), gl = 12)
    plo, phi = promote(float(_para_range(para)[1]), float(_para_range(para)[2]))
    qlo, qhi = oftype(phi, _perp_range(perp)[1]), oftype(phi, _perp_range(perp)[2])
    dg = @something dgrad (q, u) -> _grad2(f0, q, u)
    bd = isnothing(dgrad) ? (q, u) -> _val_dwrt(x -> f0(q, x), u) : (q, u) -> (f0(q, u), dgrad(q, u)[2])
    vp, up = _cross_pivots(f0, range(qlo, qhi, probe[1]), range(plo, phi, probe[2]), rtol, rmax)
    M = inv([f0(v, u) for v in vp, u in up])
    mk(cache) = LowRankVDF(
        erase_f2(f0, phi), erase_g2(dg, phi), erase_g2(bd, phi),
        vp, up, M, (plo, phi), (qlo, qhi), cache
    )
    return mk(_lr_cache(mk(nothing), gl, rtol))
end

LowRankVDF(d::CoupledVDF; kw...) = begin
    @assert regime(d) == NonRelativistic() "Only non-relativistic regime is supported"
    LowRankVDF(d.f0; para = d.para, perp = d.perp, dgrad = d.dgrad, kw...)
end

Base.show(io::IO, d::LowRankVDF) =
    print(io, "LowRankVDF(rank=", rank(d), ", para=", d.para, ", perp=", d.perp, ")")

# Adaptive cross with full pivoting on a probe grid. Rank is capped by conditioning, not
# by `rtol` alone: cond(f₀[vp,up]) ≈ 1/rtol, so the perp factors lose ~log10(1/rtol) digits.
function _cross_pivots(f0, vs, us, rtol, rmax)
    R = [float(f0(v, u)) for v in vs, u in us]
    nrm = maximum(abs, R)
    iszero(nrm) && throw(ArgumentError("LowRankVDF: f0 vanishes on the probe grid"))
    iv = Int[]; iu = Int[]
    for _ in 1:rmax
        _, idx = findmax(abs, R)
        i, j = Tuple(idx)
        abs(R[i, j]) <= rtol * nrm && break
        push!(iv, i); push!(iu, j)
        R = R .- (R[:, j] * R[i, :]') ./ R[i, j]
    end
    isempty(iv) && throw(ArgumentError("LowRankVDF: cross approximation found no pivot"))
    return collect(float.(vs[iv])), collect(float.(us[iu]))
end

# All R perp factors (and their p⊥-derivatives) at one node; one M-contraction per node.
@inline _perp_factors(d::LowRankVDF, v) =
    (d.M' * [d.f0(v, u) for u in d.up], d.M' * [d.dgrad(v, u)[1] for u in d.up])
# Parallel factor bₛ(u)=f₀(vpₛ,u) and its p∥-derivative bₛ′(u)=∂∥f₀(vpₛ,u)
@inline _b(d::LowRankVDF, s, u) = d.f0(d.vp[s], u)
@inline _bdb(d::LowRankVDF, s, u) = d.bd(d.vp[s], u)

# EVERYTHING parallel is k-independent — wavevector enters χ only through perpendicular
# Bessel argument and the pole. Keeping the p∥ table in the cache (not in the per-k plan)
# is what makes re-planning cheap.
#
# Density and ⟨p⊥²⟩ OF THE SURROGATE.
# Both are FIXED-node Gauss–Legendre on the base panels, never adaptive: the ãₛ inherit
# the cross's 1/rtol conditioning, so an adaptive rule targeting tight tolerance on them would
# chase round-off and subdivide without bound (as would a panel proxy built from them).
function _lr_cache(d::LowRankVDF, gl, rtol)
    upan = _panels(u -> sum(s -> abs(_b(d, s, u)), 1:rank(d)), d.para..., zero(d.para[2]))
    pa = ParaTable((s, u) -> _bdb(d, s, u), rank(d), upan, gl, FixedNodeEval().tol)
    vpan = _panels(v -> sum(u -> abs(d.f0(v, u)), d.up), d.perp..., zero(d.para[2]))
    n = zero(d.para[2]); p2 = zero(n)
    for (v, w) in zip(gl_nodes(vpan, QuadGK.gauss(gl)...)...)
        av = _perp_factors(d, v)[1]
        for s in 1:rank(d)
            n += 2π * w * v * av[s] * pa.mom[1, s][1]
            p2 += 2π * w * v^3 * av[s] * pa.mom[1, s][1]
        end
    end
    # ãₛ at a few perp probes, so `trust_error` costs O(nprobe + rank) f₀ evaluations per u
    vprobe = collect(range(d.perp[1], d.perp[2], _LR_NPROBE))
    A = reduce(vcat, transpose(_perp_factors(d, v)[1]) for v in vprobe)
    f0max = maximum(v -> abs(d.f0(v, zero(n))), vprobe)
    return (; n, pperp2_mean = p2 / n, para = pa, vpan, gl, rtol, vprobe, A, f0max)
end

const _LR_NPROBE = 24

"""
    trust_error(d::LowRankVDF, u)

Error of the separable expansion at (generally complex) parallel momentum `u`, 
against the true `f₀`: `max_v |f₀(v,u) − Σₛ ãₛ(v)bₛ(u)| / max|f₀|`.

Sits at `rtol` on the real axis and grows off it as `ãₛ` were fitted on the real axis.
That growth is the horizon past which a crossed-pole residue,
and hence any root behind it, is meaningless.
"""
function trust_error(d::LowRankVDF, u)
    c = d.cache
    R = rank(d)
    uc = complex(float(u))
    bs = [_b(d, s, uc) for s in 1:R]
    err = zero(real(uc))
    @inbounds for i in eachindex(c.vprobe)
        approx = zero(eltype(bs))
        for s in 1:R
            approx += c.A[i, s] * bs[s]
        end
        err = max(err, abs(d.f0(c.vprobe[i], uc) - approx))
    end
    return err / c.f0max
end

# A root is only as trustworthy as the surrogate is AT THE ζ IT SAMPLES
const _LR_TRUST = 1.0e4
function trusted(d::LowRankVDF, s, ω, k)
    kz = k.k_para
    iszero(kz) && return true
    sign(kz) * imag(ω) / kz < 0 || return true # growing side: no residue is taken
    lo, hi = d.para
    tol = _LR_TRUST * d.cache.rtol
    nmax = nmax_bessel((k.k_perp / s.Omega)^2 * abs(d.cache.pperp2_mean) / 2)
    return all((-nmax):nmax) do n
        ζ = (ω - n * s.Omega) / kz
        lo < real(ζ) < hi || return true    # pole outside the p∥ range: no residue
        trust_error(d, ζ) <= tol
    end
end

prepare(d::LowRankVDF; kw...) = PreparedVDF(d, d.cache)

plan_contribution(c::PreparedVDF{<:LowRankVDF}, s, k; kw...) =
    plan_contribution(c.vdf, s, k; kw...)

function plan_contribution(d::LowRankVDF, s, k; kw...)
    c = d.cache
    Ω = s.Omega
    a = k.k_perp / Ω
    nmax = nmax_bessel(a^2 * abs(c.pperp2_mean) / 2)
    ns = (-nmax):nmax
    # The base panels resolve the vdf; the Bessel wavelength π/a is the only k-dependent thing
    # about the perp grid (negative Ω must still refine).
    vs, ws = gl_nodes(_refine(c.vpan, abs(a)), QuadGK.gauss(c.gl)..., typeof(float(k.k_perp)))
    fac = _perp_factors.(Ref(d), vs)
    α = stack(last, fac; dims = 1)       # ãₛ′(vⱼ)
    β = stack(first, fac; dims = 1)      # ãₛ(vⱼ)
    P∂, PF = _perp_tensors(_bessel_tensors(a, ns, vs), vs, ws, α, β)
    return TablePlan(c.para, P∂, PF, ns, Ω, k.k_para, s.Pi2 / c.n)
end

contribution(d::LowRankVDF, s, ω, k; kw...) = plan_contribution(d, s, k; kw...)(ω)
contribution(c::PreparedVDF{<:LowRankVDF}, s, ω, k; kw...) = plan_contribution(c, s, k; kw...)(ω)
