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

Drop-in for [`CoupledVDF`](@ref): χ costs `O(rank · nmax)` per ω instead of a 2-D adaptive
quadrature, 10–300× faster per ω (the win grows with `k⊥`). Coupled VDFs are typically very low
rank — a bi-kappa and a spherical shell are both ≈ 10.

`f0` must be evaluable at complex `p∥`: each `bₛ(u)=f₀(vₛ,u)` is a literal `f₀` slice, so unlike a
fitted surrogate ([`GridVDF`](@ref)), which exists only on the real axis, it can be continued.

# Accuracy

The cross is fitted on the REAL axis, so `rtol` bounds χ only where the p∥ integral stays there:
real ω, growing modes, `k∥=0`.

A CROSSED Landau pole is different — it picks up the residue `f₀(p⊥,ζ)` at complex
`ζ=(ω−nΩ)/k∥`, where the `ãₛ` are extrapolating. That error grows with `|Im ω / k∥|`, is NOT
bounded by `rtol`, and eventually makes the surrogate's det gain exact zeros that are not modes of
`f₀`. Surveys drop those via [`trusted`](@ref); [`trust_error`](@ref) measures the horizon.
"""
struct LowRankVDF{F, G, H, T, C} <: AbstractVDF
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
    plo, phi = promote(float(para[1]), float(para[2]))
    qlo, qhi = oftype(phi, _pair(perp)[1]), oftype(phi, _pair(perp)[2])
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
@inline _db(d::LowRankVDF, s, u) = _bdb(d, s, u)[2]

# Neumann far-field is truncated at (1/THETA)^NMOM ≈ 1e-12.
const _LR_THETA = 2.0
const _LR_NMOM = 40
# Plan-time quadratures must not be the accuracy bottleneck: the cross tolerance `rtol` is.
const _LR_QRTOL = 1.0e-11
# Panel discovery only needs to LOCATE structure; GL-12 per panel supplies the accuracy.
const _LR_PRTOL = 1.0e-8

# EVERYTHING parallel is k-independent — wavevector enters χ only through perpendicular
# Bessel argument and the pole. Keeping the p∥ tables here (not in the per-k plan)
# is what makes re-planning cheap; `residual` re-plans once per root, and an adaptive
# moment quadrature there would dominate a whole survey.
struct LowRankPara{T}
    U::T                            # p∥ half-width: sets the far/near split
    un::Vector{T}; uw::Vector{T}    # fixed p∥ nodes/weights (near field)
    Bv::Matrix{T}; dBv::Matrix{T}   # [l,s] bₛ, bₛ′ at those nodes
    nu::Matrix{T}; mu::Matrix{T}    # [p,s] scaled p∥ moments of bₛ, bₛ′ (far field)
    I0::Vector{T}; I1::Vector{T}; I2::Vector{T}
    J0::Vector{T}; J1::Vector{T}    # raw p∥ moments (moment recursion; k∥=0 path)
end

function LowRankPara(d::LowRankVDF, gl)
    R = rank(d)
    ulo, uhi = d.para
    T = typeof(float(uhi))
    U = max(abs(ulo), abs(uhi))
    upan = _panels(u -> sum(s -> abs(_b(d, s, u)), 1:R), ulo, uhi, zero(U))
    xp, wp = QuadGK.gauss(gl)
    un = T[]; uw = T[]
    for p in 1:(length(upan) - 1)
        mid, half = (upan[p] + upan[p + 1]) / 2, (upan[p + 1] - upan[p]) / 2
        append!(un, mid .+ half .* xp); append!(uw, half .* wp)
    end
    L = length(un)
    Bv = [_b(d, s, un[l]) for l in 1:L, s in 1:R]
    dBv = [_db(d, s, un[l]) for l in 1:L, s in 1:R]
    # Far-field moments νₚ = (1/U)∫(u/U)ᵖ bₛ du need their OWN quadrature: the weight
    # (u/U)^40 is edge-peaked, so panels chosen for bₛ alone under-resolve them badly.
    nu = zeros(T, _LR_NMOM + 1, R); mu = zeros(T, _LR_NMOM + 1, R)
    for s in 1:R
        nu[:, s] .= _scaled_moments(u -> _b(d, s, u), ulo, uhi, U)
        mu[:, s] .= _scaled_moments(u -> _db(d, s, u), ulo, uhi, U)
    end
    # ∫uᵖ b du = U^{p+1}·νₚ — the moment recursion's raw moments come free.
    return LowRankPara(
        U, un, uw, Bv, dBv, nu, mu,
        U .* nu[1, :], U^2 .* nu[2, :], U^3 .* nu[3, :], U .* mu[1, :], U^2 .* mu[2, :]
    )
end

_scaled_moments(b, lo, hi, U) = QuadGK.quadgk(lo, hi; rtol = _LR_QRTOL) do u
    t = u / U
    bu = b(u) / U
    SVector(ntuple(p -> bu * t^(p - 1), Val(_LR_NMOM + 1)))
end[1]

# Density and ⟨p⊥²⟩ OF THE SURROGATE.
# Both are FIXED-node Gauss–Legendre on the base panels, never adaptive: the ãₛ inherit
# the cross's 1/rtol conditioning, so an adaptive rule targeting _LR_QRTOL on them would
# chase round-off and subdivide without bound (as would a panel proxy built from them).
function _lr_cache(d::LowRankVDF, gl, rtol)
    pa = LowRankPara(d, gl)
    vpan = _panels(v -> sum(u -> abs(d.f0(v, u)), d.up), d.perp..., zero(pa.U))
    xg, wg = QuadGK.gauss(gl)
    n = zero(pa.U); p2 = zero(n)
    for p in 1:(length(vpan) - 1)
        mid, half = (vpan[p] + vpan[p + 1]) / 2, (vpan[p + 1] - vpan[p]) / 2
        for q in eachindex(xg)
            v = mid + half * xg[q]; w = half * wg[q]
            av = _perp_factors(d, v)[1]
            for s in 1:rank(d)
                n += 2π * w * v * av[s] * pa.I0[s]
                p2 += 2π * w * v^3 * av[s] * pa.I0[s]
            end
        end
    end
    # ãₛ at a few perp probes, so `trust_error` costs O(nprobe + rank) f₀ evaluations per u
    vprobe = collect(range(d.perp[1], d.perp[2], _LR_NPROBE))
    A = reduce(vcat, transpose(_perp_factors(d, v)[1]) for v in vprobe)
    f0max = maximum(v -> abs(d.f0(v, zero(pa.U))), vprobe)
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
    kz = para(k)
    iszero(kz) && return true
    sign(kz) * imag(ω) / kz < 0 || return true # growing side: no residue is taken
    lo, hi = d.para
    tol = _LR_TRUST * d.cache.rtol
    nmax = nmax_bessel((perp(k) / s.Omega)^2 * abs(d.cache.pperp2_mean) / 2)
    return all((-nmax):nmax) do n
        ζ = (ω - n * s.Omega) / kz
        lo < real(ζ) < hi || return true    # pole outside the p∥ range: no residue
        trust_error(d, ζ) <= tol
    end
end

prepare(d::LowRankVDF, args...; kw...) = PreparedVDF(d, d.cache)

struct LowRankPlan{T, D, P}
    vdf::D
    pa::P                                 # k-independent p∥ tables
    Ω::T; kz::T; Pi2::T
    ns::UnitRange{Int}
    Pd::Matrix{SMatrix{3, 3, T, 9}}       # [n,s] 2π ∫ ãₛ′(v)·Kₙ(v) dv
    Pf::Matrix{SMatrix{3, 3, T, 9}}       # [n,s] 2π ∫ v·ãₛ(v)·Kₙ(v) dv
    invn::T
end

isexact(::LowRankPlan) = false

plan_contribution(c::PreparedVDF{<:LowRankVDF}, s, k; kw...) =
    plan_contribution(c.vdf, s, k; kw...)

function plan_contribution(d::LowRankVDF, s, k; kw...)
    c = d.cache
    gl = c.gl
    Ω, kz = s.Omega, para(k)
    a = perp(k) / Ω
    R = rank(d)
    nmax = nmax_bessel(a^2 * abs(c.pperp2_mean) / 2)
    ns = (-nmax):nmax

    # The base panels resolve the vdf; refine to ≥2 panels per Bessel wavelength π/a, which
    # is the only k-dependent thing about the perp grid (a quadrature rule sized for f₀
    # alone cannot see that oscillation).
    vpan = _refine(c.vpan, abs(a))   # Bessel wavelength π/|a|; negative Ω must still refine
    xg, wg = QuadGK.gauss(gl)
    T = typeof(float(perp(k)))
    Pd = zeros(SMatrix{3, 3, T, 9}, length(ns), R)
    Pf = zeros(SMatrix{3, 3, T, 9}, length(ns), R)
    Jv = Vector{T}(undef, nmax + 2)
    for p in 1:(length(vpan) - 1)
        mid, half = (vpan[p] + vpan[p + 1]) / 2, (vpan[p + 1] - vpan[p]) / 2
        for q in eachindex(xg)
            v = mid + half * xg[q]; w = half * wg[q]
            av, adv = _perp_factors(d, v)
            besselj_ladder!(Jv, nmax + 1, a * v)
            for (i, n) in enumerate(ns)
                Jm, Jn, Jp = _jladder(Jv, n - 1), _jladder(Jv, n), _jladder(Jv, n + 1)
                b1, b2 = v * (Jm + Jp) / 2, v * (Jm - Jp) / 2
                K = _symmat(b1^2, b1 * b2, b1 * Jn, b2^2, b2 * Jn, Jn^2)
                for r in 1:R
                    Pd[i, r] += (2π * w * adv[r]) * K
                    Pf[i, r] += (2π * w * v * av[r]) * K
                end
            end
        end
    end
    return LowRankPlan(d, c.para, Ω, kz, s.Pi2, ns, Pd, Pf, 1 / c.n)
end

# Adaptive panel edges for `f` on [lo,hi]; `a>0` additionally caps panels at the Bessel
# half-wavelength π/(2a), which adaptive quadrature on a non-oscillatory proxy cannot see.
function _panels(f, lo, hi, a, rtol = _LR_PRTOL)
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

# Landau integrals A₀ₛ=∫bₛ/(u−ζ)du, B₀ₛ=∫bₛ′/(u−ζ)du for ALL ranks at one ζ, σ=sign(k∥).
# Batched over s because ζ is shared by the whole rank index: the node kernel wₗ=uwₗ/(uₗ−ζ)
# — the only division in the ω evaluation — is built once and contracted with every column.
# FAR speedup (|ζ|>θU): Neumann series in 1/ζ — no logs, no cancellation, ω-independent moments;
#   When k∥·p∥max/Ω ≪ 1 (quasi-perpendicular propagation), all but a couple of harmonics take this branch.
# NEAR: pole-subtracted fixed-node sum + the analytic log, with bₛ(ζ),
#   bₛ′(ζ) from the TRUE f₀ slice ⇒ exact Landau residue.
function _lr_cauchy!(A, B, w, d, pa::LowRankPara, ζ, σ)
    lo, hi = d.para
    if abs(ζ) > _LR_THETA * pa.U
        invξ = pa.U / ζ
        crossed = σ * imag(ζ) < 0 && lo < real(ζ) < hi
        @inbounds for s in eachindex(A)
            a = zero(ζ); b = zero(ζ)
            for p in (_LR_NMOM + 1):-1:1
                a = (a + pa.nu[p, s]) * invξ
                b = (b + pa.mu[p, s]) * invξ
            end
            A[s] = -a; B[s] = -b
            if crossed
                bζ, dbζ = _bdb(d, s, ζ)
                A[s] += σ * _2πim * bζ
                B[s] += σ * _2πim * dbζ
            end
        end
        return
    end
    l0 = 0                      # node coinciding with a real ζ: removable, restored per s
    W = zero(ζ)
    @inbounds for l in eachindex(w)
        δ = pa.un[l] - ζ
        iszero(δ) && (l0 = l)
        w[l] = iszero(δ) ? zero(ζ) : pa.uw[l] * inv(δ)
        W += w[l]
    end
    w0 = iszero(l0) ? zero(eltype(pa.uw)) : pa.uw[l0]
    # The peel gate only needs the magnitude of φ on the real axis beside ζ, so read it off the
    # node table (bracketing the clamped ζ)
    j = searchsortedfirst(pa.un, clamp(real(ζ), lo, hi))
    j1, j2 = clamp(j - 1, 1, length(pa.un)), clamp(j, 1, length(pa.un))
    @inbounds for s in eachindex(A)
        bζ, dbζ = _bdb(d, s, ζ)
        # b″(ζ) is the removable value of the bₛ′ integral at that coincidence;
        # Note only real ζ can land on a fixed node
        d2 = iszero(l0) ? dbζ : _d2slice(d, s, real(ζ))
        sb = max(abs(pa.Bv[j1, s]), abs(pa.Bv[j2, s]))
        sdb = max(abs(pa.dBv[j1, s]), abs(pa.dBv[j2, s]))
        A[s] = _cauchy_near(w, W, w0, pa.Bv, s, bζ, dbζ, ζ, lo, hi, σ, sb)
        B[s] = _cauchy_near(w, W, w0, pa.dBv, s, dbζ, d2, ζ, lo, hi, σ, sdb)
    end
    return
end

# One near-field Landau-Cauchy transform ∫φ/(u−ζ)du from the shared kernel: Σₗwₗφₗ, minus
# φ(ζ)·Σₗwₗ when the pole is peeled, plus the analytic log.
@inline function _cauchy_near(w, W, w0, φl, s, φζ, dφζ, ζ, lo, hi, σ, scale)
    S = zero(ζ)
    @inbounds @simd for l in eachindex(w)
        S += w[l] * φl[l, s]
    end
    peeled = _peel(φζ, scale)
    peeled && (S += w0 * dφζ - φζ * W)
    return S + φζ * _lpole_term(ζ, lo, hi, σ, peeled)
end

# b″(ζ) by central difference on the first-derivative slice.
@inline _d2slice(d, s, x) = let h = 1.0e-4 * (1 + abs(x))
    (_db(d, s, x + h) - _db(d, s, x - h)) / (2h)
end

# uʲ/(u−ζ) = ζʲ/(u−ζ) + (poly) ⇒ only TWO Cauchy transforms per (n,s); the u¹,u² moments
# follow from the ω-independent raw moments.
function (pl::LowRankPlan)(ω)
    ωc = complex(float(ω))
    pa = pl.pa
    acc = zero(SVector{6, ComplexF64})
    kz = pl.kz
    z0 = iszero(kz)
    σ = z0 ? one(kz) : sign(kz)
    ik = z0 ? zero(kz) : -1 / kz
    T = typeof(ωc / oneunit(kz))
    A = similar(pa.I0, T); B = similar(pa.I0, T)   # per-rank Landau integrals at one ζ
    w = similar(pa.un, T)                          # shared near-field node kernel
    @inbounds for (i, n) in enumerate(pl.ns)
        nΩ = n * pl.Ω
        ζ = z0 ? ωc : (ωc - nΩ) / kz
        z0 || _lr_cauchy!(A, B, w, pl.vdf, pa, ζ, σ)
        for s in eachindex(pa.I0)
            M = if z0
                invΔ = 1 / (ωc - nΩ)
                (pa.I0[s] * invΔ, pa.I1[s] * invΔ, pa.I2[s] * invΔ, pa.J0[s] * invΔ, pa.J1[s] * invΔ)
            else
                A0, B0 = A[s], B[s]
                A1 = pa.I0[s] + ζ * A0
                A2 = pa.I1[s] + ζ * A1
                B1 = pa.J0[s] + ζ * B0
                (ik * A0, ik * A1, ik * A2, ik * B0, ik * B1)
            end
            acc += _chi_mblock(M, pl.Pd[i, s], pl.Pf[i, s], ωc, kz, nΩ)
        end
    end
    return (pl.Pi2 * pl.invn / ωc^2) * _antisymmat(acc)
end

contribution(d::LowRankVDF, s, ω, k; kw...) = plan_contribution(d, s, k; kw...)(ω)
contribution(c::PreparedVDF{<:LowRankVDF}, s, ω, k; kw...) = plan_contribution(c, s, k; kw...)(ω)
