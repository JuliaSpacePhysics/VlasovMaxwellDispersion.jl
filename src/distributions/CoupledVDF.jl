"""
    CoupledVDF(f0; para=(-Inf,Inf), perp=Inf, dgrad=nothing, n=nothing, regime=NonRelativistic())

**Most general** gyrotropic VDF: an arbitrary analytic `f0(p⊥,p∥)`.

`f0` must be evaluable at complex argument (continued onto the Landau contour) and,
with the default infinite bounds, safe at arbitrarily large real argument (return 0,
don't overflow).

It is stored raw; χ is linear in `f₀`, so `contribution` scales it by `1/n`.

`para`/`perp` integration bounds are optional: the default integrates p∥ over ℝ and
p⊥ over (0,∞) via a sinh-mapped trapezoid with an exact cotangent pole correction
(see `_coupled_perp_sinc`) — bounds-free, exponentially convergent, and typically
faster than a hand-tuned box. Pass finite `para=(lo,hi)`, `perp=(lo,hi)` to force the
adaptive finite-box path (required for `regime=Relativistic()` and
`closure=Newberger()`, and the right choice for piecewise-smooth `f0` such as spline
fits, where the trapezoid loses its spectral accuracy). Bounds must be all finite or
all infinite.

`dgrad(p⊥,p∥) -> (∂⊥f0, ∂∥f0)` supplies the gradient and default to autodiff.

`n ≡ ∫d³p f₀`; default `nothing` autocomputes `n` by quadrature.

Prefer [`SeparableVDF`] when `f0(p⊥,p∥)=f⊥(p⊥)f∥(p∥)`.
"""
struct CoupledVDF{F, Dg, T, R <: Regime, B} <: AbstractVDF
    f0::F
    dgrad::Dg
    para::Tuple{T, T}
    perp::Tuple{T, T}
    n::T           # ∫d³p f₀
    pperp2::T      # ⟨p⊥²⟩: Bessel harmonic window (nmax)
    upar::T        # ⟨p∥⟩: sinc-map center
    ppar2::T       # ⟨(p∥−⟨p∥⟩)²⟩: sinc-map scale
    regime::R
    # B (finite bounds?) is a TYPE parameter so the box and sinc harmonic paths are
    # separate specializations — a runtime branch would compile both for every f0.
    function CoupledVDF(
            f0::F, dgrad::Dg, para::Tuple{T, T}, perp::Tuple{T, T}, n::T,
            pperp2::T, upar::T, ppar2::T, regime::R
        ) where {F, Dg, T, R <: Regime}
        B = isfinite(para[1]) & isfinite(para[2]) & isfinite(perp[2])
        return new{F, Dg, T, R, B}(f0, dgrad, para, perp, n, pperp2, upar, ppar2, regime)
    end
end

const FiniteCoupledVDF = CoupledVDF{<:Any, <:Any, <:Any, <:Any, true}
const InfiniteCoupledVDF = CoupledVDF{<:Any, <:Any, <:Any, <:Any, false}

regime(d::CoupledVDF) = d.regime

@inline _pair(x::Tuple) = x
@inline _pair(x) = (zero(x), x)

function CoupledVDF(
        f0; para = (-Inf, Inf), perp = Inf, dgrad = nothing, n = nothing,
        regime = NonRelativistic()
    )
    plo, phi = promote(float(para[1]), float(para[2]))
    qlo, qhi = oftype(phi, _pair(perp)[1]), oftype(phi, _pair(perp)[2])
    finite = isfinite(plo) & isfinite(phi) & isfinite(qhi)
    finite || (isinf(plo) & isinf(phi) & isinf(qhi)) ||
        throw(ArgumentError("CoupledVDF bounds must be all finite (box) or all infinite (bounds-free default)"))
    regime isa Relativistic && !finite &&
        throw(ArgumentError("relativistic CoupledVDF needs finite para/perp bounds (fixed-order box quadrature)"))
    n = @something n 2π * QuadGK.quadgk(
        q -> q * QuadGK.quadgk(u -> f0(q, u), plo, phi; rtol = 1.0e-9)[1],
        qlo, qhi; rtol = 1.0e-9
    )[1]
    # f0-only moments, hoisted so χ calls skip the pre-integrals. ⟨p∥⟩/⟨p∥²⟩ must exist
    # for the default infinite bounds (heavy tails need κ>3/2-type decay). MUST stay one
    # fused vector quadrature: ⟨p∥⟩ is exactly zero for symmetric f0, and a scalar quadgk
    # can never meet a relative tolerance on a zero integral (refines to maxevals in
    # every nested inner call — constructor hang); the vector norm lends the zero
    # component the nonzero components' scale.
    mom = 2π * QuadGK.quadgk(qlo, qhi; rtol = 1.0e-3) do q
        inner = QuadGK.quadgk(plo, phi; rtol = 1.0e-3) do u
            fv = f0(q, u)
            SVector(fv, u * fv, u^2 * fv)
        end[1]
        SVector(q^3 * inner[1], q * inner[2], q * inner[3])
    end[1] / n
    dg = isnothing(dgrad) ? ((q, u) -> _grad2(f0, q, u)) : dgrad
    return CoupledVDF(
        f0, dg, (plo, phi), (qlo, qhi), oftype(phi, n),
        oftype(phi, abs(mom[1])), oftype(phi, mom[2]), oftype(phi, abs(mom[3] - mom[2]^2)), regime
    )
end

function contribution(d::CoupledVDF, s, ω, k; closure = HarmonicSum(), kwargs...)
    return _coupled_contribution(closure, regime(d), d, s, complex(float(ω)), k; kwargs...) / d.n
end

function _coupled_contribution(::HarmonicSum, ::NonRelativistic, d::FiniteCoupledVDF, s, ω, k; norm = NORM, rtol = 1.0e-6)
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    a = kperp / Ω
    L, U = d.para
    nmax = nmax_bessel(a^2 * d.pperp2 / 2)
    ns = (-nmax):nmax
    ζs = [(ω - n * Ω) / kz for n in ns]
    X = QuadGK.quadgk(d.perp...; rtol, norm) do v
        _coupled_perp(v, ns, ζs, d, ω, Ω, kz, a, L, U; norm, rtol)
    end[1]
    return (s.Pi2 / ω^2) * _antisymmat(X)
end

# Bounds-free path: sinh-mapped trapezoid + cotangent pole correction — see
# experiments/infinite-bounds/README.md for derivation and benchmarks.
function _coupled_contribution(::HarmonicSum, ::NonRelativistic, d::InfiniteCoupledVDF, s, ω, k; norm = NORM, rtol = 1.0e-6, h = 0.2)
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    a = kperp / Ω
    nmax = nmax_bessel(a^2 * d.pperp2 / 2)
    ns = (-nmax):nmax
    ζs = [(ω - n * Ω) / kz for n in ns]
    uc, S = d.upar, sqrt(d.ppar2)
    S > 0 || throw(ArgumentError("CoupledVDF sinc path: ⟨p∥²⟩ vanished — f0 has no parallel width"))
    X = QuadGK.quadgk(d.perp...; rtol, norm) do v
        _coupled_perp_sinc(v, ns, ζs, d, ω, Ω, kz, a, uc, S, h; norm)
    end[1]
    return (s.Pi2 / ω^2) * _antisymmat(X)
end

# Landau-continued parallel Cauchy integrals for the WHOLE harmonic ladder at one perp
# node, by Kress-type product quadrature: with u = ψ(t) = uc + S·sinh(t) and uniform
# nodes t_j = j·h (alignment to j·h is REQUIRED — an offset grid breaks the cot phase),
#
#   ∫_ℝ g(u)/(u−ζ) du = h·Σ_j G(t_j)/(ψ(t_j)−ζ) + π·g(ζ)·(cot(π·t*/h) + i),
#   G = g(ψ)ψ′,  t* = asinh((ζ−uc)/S),
#
# from residue calculus on (π/h)cot(πw/h)·G/(ψ−ζ): the residue at t* is g(ζ) exactly
# (Jacobian cancels), and the single +iπ merges the Im ζ ≷ 0 cases WITH the Landau 2πi
# — one analytic formula, no crossed bookkeeping, exact-real ζ included. Error
# ~e^{−2πd/h} (d = analyticity strip of G): h=0.25→1e-7, 0.2→1e-9, 0.15→1e-13.
# All harmonics share the fixed samples; nodes where f0 ≈ 0 skip the gradient evals.
function _coupled_perp_sinc(v, ns, ζs, d::CoupledVDF, ω, Ω, kz, a, uc, S, h; T = 7.0, ftol = 1.0e-14, kw...)
    g5(u) = begin
        q, p = d.dgrad(v, u)
        SVector(q, u * q, u^2 * q, p, u * p)
    end
    invkz = -1 / kz
    nb = length(ns)
    jmax = floor(Int, T / h)
    return @no_escape begin
        b2s = @alloc(SVector{6, typeof(a * v)}, nb)
        _perp_Bessel_bilinears!(b2s, a, v)
        Is = @alloc(SVector{5, eltype(ζs)}, nb)
        @inbounds for i in 1:nb
            Is[i] = zero(SVector{5, eltype(ζs)})
        end
        thr = ftol * (abs(d.f0(v, uc)) + abs(d.f0(v, uc + S)) + abs(d.f0(v, uc - S)))
        for j in (-jmax):jmax
            t = j * h
            u = uc + S * sinh(t)
            abs(d.f0(v, u)) < thr && continue
            w = h * S * cosh(t)
            g = g5(u)
            @inbounds for i in 1:nb
                Is[i] += (w / (u - ζs[i])) * g
            end
        end
        acc = zero(AType)
        @inbounds for i in 1:nb
            gz = g5(ζs[i])
            corr = all(isfinite, gz) ? (π * _cot_i(π * asinh((ζs[i] - uc) / S) / h)) * gz :
                (imag(ζs[i]) < 0 ? gz .* (2π * im) : zero(gz))
            acc += _In_block(Is[i] + corr, invkz, b2s[i], v, ω, kz, ns[i] * Ω)
        end
        acc
    end
end

# cot(w) + i, saturated: cot → ∓i as Im w → ±∞ but overflows past |Im w| ≈ 700;
# beyond 20 the residual is < 4e-18, so return the limit (0 or 2i) exactly.
@inline _cot_i(w) = abs(imag(w)) > 20 ? complex(0.0, imag(w) > 0 ? 0.0 : 2.0) : cot(w) + im

# Relativistic path, sliced in (p⊥,p∥) — docs/relativistic.md.
# Resonance D(p∥) = ωγ − k∥p∥ − nΩ₀ with γ=√(1+p⊥²+p∥²) rationalizes,
#   D·D̃ = A(p∥−p₊)(p∥−p₋),  D̃ = ωγ + k∥p∥ + nΩ₀,  A = ω²−k∥²,
# into two explicit simple poles; the squaring ghost (zero of D̃) carries a null
# residue automatically. Poles cross the real p∥ axis ONLY at Im ω = 0.
# Endpoints |p∥|=P sit where f₀≈0, so no endpoint (rim-type) corrections arise.
# f₀ must be evaluable at complex p∥ (poles sit off-axis for complex ω).
# Validated vs Maxwell–Jüttner (Swanson) to ~1e-5 down to Im ω = −0.15 at μ=2.
function _coupled_contribution(::HarmonicSum, ::Relativistic, d::CoupledVDF, s, ω, k; rtol = 1.0e-6)
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    if imag(ω) < 0 && real(ω)^2 > kz^2
        @warn "damped superluminal ω (|Re ω| > |k∥|): the (p⊥,p∥) integral is not the analytic continuation there (apex branch cut, docs/relativistic.md); evaluate at Im ω ≥ 0 and continue externally" maxlog = 1
    end
    a = kperp / Ω
    plo, phi = d.para
    qhi = d.perp[2]
    γmax = sqrt(1 + max(phi^2, plo^2) + qhi^2)
    nmax = nmax_bessel(a^2 * qhi^2 / 2)
    f = n -> _coupled_harmonic_rel(n, d, ω, Ω, kz, a)
    X_T = 2π * converge(f; nmax, rtol)
    X = _antisymmat(X_T) .+ _ee33(_bernstein_rel(d, γmax))
    return (s.Pi2 / ω^2) * X
end

# Relativistic non-resonant e∥e∥ term without prefactor
function _bernstein_rel(d, γmax; GLγ = _GLγ, GLp = _GLp)
    gn, gw = GLγ
    pn, pw = GLp
    acc = zero(ComplexF64)
    for ig in eachindex(gn)
        q = (gn[ig] + 1) / 2
        γ = 1 + (γmax - 1) * q^2
        wγ = gw[ig] * (γmax - 1) * q
        umax = sqrt(γ^2 - 1)
        inner = zero(ComplexF64)
        for ip in eachindex(pn)
            θ = pn[ip] * (π / 2)
            u, w = umax .* sincos(θ)
            dpe, dpa = d.dgrad(w, u)
            inner += pw[ip] * (π / 2) * ComplexF64(w * u * dpa - u^2 * dpe)
        end
        acc += wγ * inner
    end
    return 2π * acc
end

@inline _ee33(x) = @SMatrix ComplexF64[0 0 0; 0 0 0; 0 0 x]


# Fixed Gauss–Legendre orders for the edge-mapped relativistic path (outer γ→q, inner p∥→θ).
# Very sharp/multi-scale f₀ may need higher orders — bump these.
const _GLγ = QuadGK.gauss(24)
const _GLp = QuadGK.gauss(32)

# Covariant momentum numerator 𝒰 = ω∂_γf+k∥∂_uf at (γ,p∥) with w=p⊥, rewritten via
# ∂_γ|_u=(γ/w)∂_⊥, ∂_u|_γ=∂_∥−(u/w)∂_⊥ ⇒ 𝒰 = k∥∂_∥f + (ωγ−k∥u)/w · ∂_⊥f.
@inline function _U_cov(d, u, w, γ, ω, kz)
    dpe, dpa = d.dgrad(w, u)
    return kz * dpa + dpe * (ω * γ - kz * u) / w
end


const AType = SVector{6, ComplexF64}


# Poles farther than this from the real p∥ segment leave the integrand smooth enough
# for plain adaptive quadrature; nearer (or Landau-crossed) poles are peeled. Also keeps
# the γ=0 artifact roots of D·D̃ (at p∥=±i√(1+p⊥²), |Im|≥1) out of the peeled set.
const _PQ_NEAR = 1.5

# Partial fractions of the rationalized resonance (docs/relativistic.md) at fixed
# m⊥² = 1+p⊥²: D_n·D̃_n = A·u² + B·u + C, so
#   1/D_n = D̃_n·[c₁/(u−p₁) + c₂/(u−p₂)],  c₁₂ = ∓1/√(B²−4AC).
# Vieta gives the second root without cancellation; A→0 sends p₁→∞, marked non-finite
# (its term is O(A); callers drop non-finite poles).
@inline function _Dn_poles(ω, kz, nΩ, m2)
    A = ω^2 - kz^2
    B = -2 * kz * nΩ
    C = ω^2 * m2 - nΩ^2
    sq = sqrt(B^2 - 4 * A * C)
    abs2(B + sq) < abs2(B - sq) && (sq = -sq)
    return (_home_side((-B - sq) / (2A), ω, kz, m2), -1 / sq),
        (_home_side(2 * C / (-B - sq), ω, kz, m2), 1 / sq)
end

# Exactly-real ω leaves a real pole ON the path: a signed zero nudges it to its home
# side (the Im ω→0⁺ limit, slope dp/dω = γ²/(k∥γ−ωp)) so the boundary-value log in
# `_peel_residue` lands on the correct sheet.
@inline function _home_side(p, ω, kz, m2)
    isfinite(p) || return complex(Inf)
    iszero(imag(p)) || return p
    γp = sqrt(complex(m2 + p^2))
    return complex(real(p), sign(real(γp^2 / (kz * γp - ω * p))) * 0.0)
end

# Residue r = c·W(p) of a peeled pole and its analytic across-box term
# r·[log((hi−p)/(lo−p)) + 2πi if Landau-crossed]; (0, 0) when not peeled.
# Peel when Landau-crossed (Im ω<0 dragged the pole below the axis inside the box —
# the +2πi is the continuation residue) or within _PQ_NEAR of the segment (Plemelj
# subtraction for quadrature health). The squaring ghost peels harmlessly (W(p)=0 ⇒
# r=0); γ-artifact roots (γ(p)=0 ⇒ W=∞) are left unpeeled.
@inline function _peel_residue(p, c, W, γof, ν, lo, hi)
    zz = (zero(AType), zero(AType))
    isfinite(p) || return zz
    crossed = ν < 0 && imag(p) < 0 && lo < real(p) < hi
    near = abs(imag(p)) < _PQ_NEAR && lo - _PQ_NEAR < real(p) < hi + _PQ_NEAR
    (crossed || near) || return zz
    r = c .* W(p, γof(p))
    all(isfinite, r) || return zz
    return r, r .* (log((hi - p) / (lo - p)) + (crossed ? 2π * im : 0))
end

# One harmonic of the (p⊥,p∥) box integral: outer Gauss–Legendre in p⊥; per slice
#   ∫ σ𝓣_n/D_n du = ∫ Σᵢ cᵢ·W/(u−pᵢ) du,  W = σ·D̃_n·𝓣_n,  σ = 𝒰·p⊥/γ,
# with peeled poles kept as single fractions (c·W(u)−r)/(u−p) — the split form
# σ𝓣/D − r/(u−p) carries 1/(u−p)² rounding noise near the pole.
function _coupled_harmonic_rel(n, d, ω, Ω, kz, a; GLq = _GLγ, GLp = _GLp)
    plo, phi = d.para
    qlo, qhi = d.perp
    nΩ = n * Ω
    ν = imag(ω)
    qn, qw = GLq
    un, uw = GLp
    qmid, qhalf = (qlo + qhi) / 2, (qhi - qlo) / 2
    umid, uhalf = (plo + phi) / 2, (phi - plo) / 2
    total = zero(AType)
    for iq in eachindex(qn)
        q = qmid + qhalf * qn[iq]
        m2 = 1 + q^2
        z = a * q
        γof(u) = sqrt(complex(m2 + u^2))
        σof(u, γ) = begin
            dpe, dpa = d.dgrad(q, u)
            (kz * dpa + dpe * (ω * γ - kz * u) / q) * (q / γ)
        end
        Wof = (u, γ) -> (σof(u, γ) * (ω * γ + kz * u + nΩ)) .* _T_n_bare(n, z, u, q)
        (p1, c1), (p2, c2) = _Dn_poles(ω, kz, nΩ, m2)
        r1, lg1 = _peel_residue(p1, c1, Wof, γof, ν, plo, phi)
        r2, lg2 = _peel_residue(p2, c2, Wof, γof, ν, plo, phi)
        reg = zero(AType)
        for iu in eachindex(un)
            u = umid + uhalf * un[iu]
            γ = γof(u)
            acc = zero(AType)
            if isfinite(p1) || isfinite(p2)
                Wu = Wof(u, γ)
                isfinite(p1) && (acc = acc .+ (c1 .* Wu .- r1) ./ (u - p1))
                isfinite(p2) && (acc = acc .+ (c2 .* Wu .- r2) ./ (u - p2))
            else                       # A=B=0 (ω=±k∥, n=0): quadratic degenerate, no poles
                acc = (σof(u, γ) .* _T_n_bare(n, z, u, q)) ./ (ω * γ - kz * u - nΩ)
            end
            reg = reg .+ (uhalf * uw[iu]) .* acc
        end
        total = total .+ (qhalf * qw[iq]) .* (reg .+ lg1 .+ lg2)
    end
    return total
end

# I(p⊥) for the WHOLE harmonic sum at one perp node
function _coupled_perp(v, ns, ζs, d::CoupledVDF, ω, Ω, kz, a, L, U; kw...)
    g5(u) = begin
        q, p = d.dgrad(v, u)
        SVector(q, u * q, u^2 * q, p, u * p)
    end
    invkz = -1 / kz
    nb = length(ns)
    gscale = maximum(ζ -> _relsize(g5(clamp(real(ζ), L, U))), ζs)
    return @no_escape begin
        gζs = @alloc(SVector{5, eltype(ζs)}, nb)
        near = @alloc(Bool, nb)
        @inbounds for i in 1:nb
            gζs[i] = g5(ζs[i])
            near[i] = _subtract_safe(gζs[i], gscale)
        end
        b2s = @alloc(SVector{6, typeof(a * v)}, nb)
        _perp_Bessel_bilinears!(b2s, a, v)
        reg = QuadGK.quadgk(L, U; kw...) do u
            g = g5(u)
            acc = zero(AType)
            @inbounds for i in eachindex(ns)
                c = invkz / (u - ζs[i])
                acc += _In_block(near[i] ? g - gζs[i] : g, c, b2s[i], v, ω, kz, ns[i] * Ω)
            end
            acc
        end[1]
        # analytic pole term, constant in u
        logacc = zero(AType)
        @inbounds for i in eachindex(ns)
            logacc += _In_block(_pole_corr(near[i], gζs[i], ζs[i], L, U), invkz, b2s[i], v, ω, kz, ns[i] * Ω)
        end
        reg + logacc
    end
end

include("qin.jl")
