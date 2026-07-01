"""
    CoupledVDF(f0; para=(lo,hi), perp=(lo,hi), dgrad=nothing, regime=NonRelativistic())

**Most general** gyrotropic VDF: an arbitrary analytic `f0(p⊥,p∥)`.

`f0` must be evaluable at complex argument (continued onto the Landau contour).

And `para`/`perp` are `(lower, upper)` integration ranges.

`dgrad(p⊥,p∥) -> (∂⊥f0, ∂∥f0)` supplies the gradient and default to autodiff.

`regime` type picks the coordinate system:
- `NonRelativistic` (default) → (p⊥,p∥)
- `Relativistic` → (γ,p∥)

Prefer [`SeparableVDF`] when `f0(p⊥,p∥)=f⊥(p⊥)f∥(p∥)`.
"""
struct CoupledVDF{F, Dg, T, R <: Regime} <: AbstractVDF
    f0::F
    dgrad::Dg
    para::Tuple{T, T}
    perp::Tuple{T, T}
    regime::R
end

regime(d::CoupledVDF) = d.regime

@inline _pair(x::Tuple) = x
@inline _pair(x) = (zero(x), x)

function CoupledVDF(
        f0; para, perp, dgrad = nothing, normalize = true,
        regime = NonRelativistic()
    )
    plo, phi = promote(para[1], para[2])
    qlo, qhi = oftype(phi, _pair(perp)[1]), oftype(phi, _pair(perp)[2])
    n = normalize ?
        2π * QuadGK.quadgk(
            q -> q * QuadGK.quadgk(u -> f0(q, u), plo, phi; rtol = 1.0e-9)[1],
            qlo, qhi; rtol = 1.0e-9
        )[1] : one(plo)
    fn = (q, u) -> f0(q, u) / n
    dg = isnothing(dgrad) ? ((q, u) -> _grad2(fn, q, u)) : ((q, u) -> dgrad(q, u) ./ n)
    return CoupledVDF(fn, dg, (plo, phi), (qlo, qhi), regime)
end

function contribution(d::CoupledVDF, s, ω, k; closure = HarmonicSum())
    return _coupled_contribution(closure, regime(d), d, s, complex(float(ω)), k)
end

function _coupled_contribution(::HarmonicSum, ::NonRelativistic, d::CoupledVDF, s, ω, k; norm = NORM, rtol = 1.0e-7)
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    a = kperp / Ω
    L, U = d.para
    p⊥²_mean = 2π * QuadGK.quadgk(
        v -> v^3 * QuadGK.quadgk(u -> d.f0(v, u), L, U; rtol)[1],
        d.perp...; rtol
    )[1]
    nmax = nmax_bessel(a^2 * abs(p⊥²_mean) / 2)
    ns = (-nmax):nmax
    ζs = [(ω - n * Ω) / kz for n in ns]
    X = QuadGK.quadgk(d.perp...; rtol, norm) do v
        _coupled_perp(v, ns, ζs, d, ω, Ω, kz, a, L, U; norm, rtol)
    end[1]
    return (s.Pi2 / ω^2) * _antisymmat(X)
end

# Relativistic (γ,p∥) momentum-space path. Momentum distribution f₀ must be
# evaluable at complex p⊥ (the pole pushes p⊥ off-axis).
# Validated vs Maxwell–Jüttner (Swanson) to ~1e-6 and → bi-Maxwellian as μ→∞.
function _coupled_contribution(::HarmonicSum, ::Relativistic, d::CoupledVDF, s, ω, k; rtol = 1.0e-6)
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    a = kperp / Ω
    plo, phi = d.para
    qhi = d.perp[2]
    γmax = sqrt(1 + max(phi^2, plo^2) + qhi^2)
    nmax = nmax_bessel(a^2 * qhi^2 / 2)
    f = n -> _coupled_harmonic_rel(n, d, ω, Ω, kz, a, γmax)
    X_T = converge(f; nmax, rtol)
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

# 3×3 relativistic harmonic integrand 2π·𝒰·𝓣_n at (γ,p∥); bare momenta make 𝓣_n
# regular at w=0. Caller passes w=√(γ²−1−u²) (complex off the real p∥ range).
@inline _rel_integrand(u, w, γ, n, a, ω, kz, d) = (2π * _U_cov(d, u, w, γ, ω, kz)) .* _T_n_bare(n, a * w, u, w)
@inline _rel_integrand(u, γ, n, a, ω, kz, d) = _rel_integrand(u, sqrt(complex(γ^2 - 1 - u^2)), γ, n, a, ω, kz, d)

const AType = SVector{6, ComplexF64}

# One relativistic harmonic, edge-mapped (derivation §5.2.2).
# Map the disk (γ,p∥) → fixed box (q,θ)∈[0,1]×[−π/2,π/2]:
#   p∥=umax·sinθ, p⊥=umax·cosθ  — inner Jacobian p⊥ cancels the rim 1/p⊥ exactly;
#   γ=1+(γmax−1)q²              — outer Jacobian ∝q flattens the √(γ−1) floor.
# Bessel stays on the fast real path.
# Off-disk poles (this n doesn't resonate at this γ) aren't peeled — nζ=0 there, so the subtraction reduces to direct integration
function _coupled_harmonic_rel(n, d, ω, Ω, kz, a, γmax; GLγ = _GLγ, GLp = _GLp)
    gn, gw = GLγ
    pn, pw = GLp
    acc = zero(AType)
    for ig in eachindex(gn)
        q = (gn[ig] + 1) / 2
        γ = 1 + (γmax - 1) * q^2
        wγ = gw[ig] * (γmax - 1) * q             # gw·½·2(γmax−1)q
        umax = sqrt(γ^2 - 1)
        ζ = (γ * ω - n * Ω) / kz                 # single Landau pole in p∥
        inrange = -umax < real(ζ) < umax
        nζ = inrange ? _rel_integrand(ζ, γ, n, a, ω, kz, d) : zero(AType)
        inner = zero(AType)
        for ip in eachindex(pn)
            θ = pn[ip] * (π / 2)
            u, w = umax .* sincos(θ)  # p⊥=w real on the disk
            wu = pw[ip] * (π / 2) * w             # Jacobian p⊥·dθ cancels rim 1/p⊥
            inner = inner .+ wu .* ((_rel_integrand(u, w, γ, n, a, ω, kz, d) .- nζ) ./ (u - ζ))
        end
        inrange && (inner = inner .+ nζ .* _landau_logfac(ζ, -umax, umax))
        acc = acc .+ wγ .* ((-1 / kz) .* inner)
    end
    return acc
end

# I(p⊥) for the WHOLE harmonic sum at one perp node
function _coupled_perp(v, ns, ζs, d::CoupledVDF, ω, Ω, kz, a, L, U; kw...)
    g5(u) = begin
        q, p = d.dgrad(v, u)
        SVector(q, u * q, u^2 * q, p, u * p)
    end
    gζs = g5.(ζs)
    b2s = _perp_Bessel_bilinears(ns, a, v)
    invkz = -1 / kz
    reg = QuadGK.quadgk(L, U; kw...) do u
        g = g5(u)
        acc = zero(AType)
        @inbounds for i in eachindex(ns)
            c = invkz / (u - ζs[i])
            acc += _In_block(g - gζs[i], c, b2s[i], v, ω, kz, ns[i] * Ω)
        end
        acc
    end[1]
    # analytic log-ratio (+ Landau) part, constant in u
    logacc = zero(AType)
    @inbounds for i in eachindex(ns)
        logacc += _In_block(gζs[i] .* _landau_logfac(ζs[i], L, U), invkz, b2s[i], v, ω, kz, ns[i] * Ω)
    end
    return reg + logacc
end

include("qin.jl")
