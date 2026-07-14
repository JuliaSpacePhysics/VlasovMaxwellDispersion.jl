# Motivation: At fixed k the ONLY ω-dependence of expensive integral in χ comes from the pole ζₙ=(ω−nΩ)/k∥;
# the perp Bessel is ω-independent. For a coupled f₀ the inner Cauchy transform still
# depends on p⊥, so perp quadrature cannot leave the ω loop.
# A rank-R separable approximation
#     f₀(p⊥,p∥) ≈ Σₛ ãₛ(p⊥)·bₛ(p∥)
# breaks it: the perp Bessel moments become ω-independent tensors P[n,s] built once
# per k, and each ω costs only R·(2nmax+1) scalar Cauchy transforms.

"""
    LowRankVDF(f0; para, perp, dgrad=nothing, rtol=1e-8, rmax=40, probe=(200, 401))
    LowRankVDF(d::CoupledVDF; kw...)

Separable surrogate `f₀ ≈ Σₛ ãₛ(p⊥)·bₛ(p∥)` of a general coupled gyrotropic `f₀(p⊥,p∥)`,
built by adaptive cross approximation to relative tolerance `rtol` (rank capped by `rmax`,
pivots searched on a `probe = (nperp, npara)` grid).

A drop-in replacement for [`CoupledVDF`](@ref) whose susceptibility costs
`O(rank · nmax)` per ω instead of a 2-D adaptive quadrature. Typical coupled VDFs are very
low rank (a spherical shell and a bi-kappa are both ≈ 10), giving 10–300× per-`ω` speedups —
the win grows with `k⊥` (harmonic count).

Approximate by construction; `rtol` sets the susceptibility's relative accuracy.

`f0` must be evaluable at complex `p∥`: each retained parallel factor `bₛ(u)=f₀(vₛ,u)` is a
literal `f₀` slice, analytic, so its Landau residue is exact.

A fitted surrogate (spline like [`GridVDF`](@ref) or SVD vector) only
exists on the real axis and cannot be continued, giving wrong damping rates.

Note that the cross truncates on the REAL axis: structure below `rtol` there is
dropped even if it dominates off-axis after continuation, so deeply damped accuracy is NOT
bounded by `rtol` (a perturbation `f₀(1+εq²cos50u)`, `ε≪rtol`, is invisible to the cross yet
sizeable at `Im ω<0`). Validate the growth rate against [`CoupledVDF`](@ref) before trusting it.
"""
struct LowRankVDF{F, G, T, C} <: AbstractVDF
    f0::F
    dgrad::G
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
    vp, up = _cross_pivots(f0, range(qlo, qhi, probe[1]), range(plo, phi, probe[2]), rtol, rmax)
    M = inv([f0(v, u) for v in vp, u in up])
    mk(cache) = LowRankVDF(erase_f2(f0, phi), erase_g2(dg, phi), vp, up, M, (plo, phi), (qlo, qhi), cache)
    return mk(_lr_cache(mk(nothing), gl))
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
@inline _db(d::LowRankVDF, s, u) = d.dgrad(d.vp[s], u)[2]

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
function _lr_cache(d::LowRankVDF, gl)
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
    return (; n, pperp2_mean = p2 / n, para = pa, vpan, gl)
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

# Cauchy transforms A₀=∫bₛ/(u−ζ)du, B₀=∫bₛ′/(u−ζ)du on the Landau sheet, σ=sign(k∥).
# FAR (|ζ|>θU): Neumann series in 1/ζ — no logs, no cancellation, ω-independent moments;
#   a Landau-crossed pole (Re ζ∈(lo,hi) dragged to the damped side) still adds its residue
#   σ·2πi·bₛ(ζ). At quasi-perpendicular propagation k∥·p∥max/Ω ≪ 1 all but a couple of
#   harmonics take this branch — the source of the speedup.
# NEAR: pole-subtracted fixed-node sum + the analytic log, with bₛ(ζ),
#   bₛ′(ζ) from the TRUE f₀ slice ⇒ exact Landau residue.
@inline function _lr_cauchy(d, pa::LowRankPara, s, ζ, σ)
    lo, hi = d.para
    if abs(ζ) > _LR_THETA * pa.U
        invξ = pa.U / ζ
        A = zero(ζ); B = zero(ζ)
        @inbounds for p in _LR_NMOM:-1:0
            A = (A + pa.nu[p + 1, s]) * invξ
            B = (B + pa.mu[p + 1, s]) * invξ
        end
        A = -A; B = -B
        if σ * imag(ζ) < 0 && lo < real(ζ) < hi
            A += σ * _2πim * _b(d, s, ζ)
            B += σ * _2πim * _db(d, s, ζ)
        end
        return (A, B)
    end
    bζ = _b(d, s, ζ)
    dbζ = _db(d, s, ζ)
    uc = clamp(real(ζ), lo, hi)
    # b″(ζ) is the removable value of the bₛ′ integral at a real-node coincidence; only real ζ
    # can land on a fixed node, so compute it only there (nested AD would collide the HoloTag).
    d2 = iszero(imag(ζ)) ? _d2slice(d, s, real(ζ)) : dbζ
    A = _cauchy_near(pa.un, pa.uw, pa.Bv, s, bζ, dbζ, ζ, lo, hi, σ, abs(_b(d, s, uc)))
    B = _cauchy_near(pa.un, pa.uw, pa.dBv, s, dbζ, d2, ζ, lo, hi, σ, abs(_db(d, s, uc)))
    return (A, B)
end

# One near-field Landau-Cauchy transform ∫φ/(u−ζ)du by fixed-node quadrature. `φl[:,s]` samples
# φ at the nodes, `φζ`=φ(ζ), `dφζ`=φ′(ζ) — the L'Hôpital value used when a node coincides with a
# real ζ (the 0*safe_inv(0) removable singularity).
@inline function _cauchy_near(un, uw, φl, s, φζ, dφζ, ζ, lo, hi, σ, scale)
    peeled = _peel(φζ, scale)
    realζ = iszero(imag(ζ))
    A = zero(ζ)
    @inbounds for l in eachindex(un)
        δ = un[l] - ζ
        A += if !peeled
            uw[l] * φl[l, s] * safe_inv(δ)
        elseif realζ && iszero(real(δ))
            uw[l] * dφζ
        else
            uw[l] * (φl[l, s] - φζ) * safe_inv(δ)
        end
    end
    return A + φζ * _lpole_term(ζ, lo, hi, σ, peeled)
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
    @inbounds for (i, n) in enumerate(pl.ns)
        nΩ = n * pl.Ω
        ζ = z0 ? ωc : (ωc - nΩ) / kz
        for s in eachindex(pa.I0)
            M = if z0
                invΔ = 1 / (ωc - nΩ)
                (pa.I0[s] * invΔ, pa.I1[s] * invΔ, pa.I2[s] * invΔ, pa.J0[s] * invΔ, pa.J1[s] * invΔ)
            else
                A0, B0 = _lr_cauchy(pl.vdf, pa, s, ζ, σ)
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
