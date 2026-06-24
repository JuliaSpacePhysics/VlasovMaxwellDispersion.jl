"""
    CoupledVDF(f0; parlower, parupper, perpupper, dpar=nothing, dperp=nothing)

**Most general** gyrotropic VDF: an arbitrary analytic `f0(p∥,p⊥)`.

`f0` must be evaluable at complex argument (continued onto the Landau contour).

Relativistic species (`regime=Relativistic()`) integrate in `(γ,p∥)`.

Prefer [`SeparableVDF`] when `f0(p∥,p⊥)=f∥(p∥)f⊥(p⊥)`.
"""
struct CoupledVDF{F,Dp,Dq,T} <: AbstractVDF
    f0::F
    dpar::Dp        # ∂f₀/∂p∥
    dperp::Dq       # ∂f₀/∂p⊥
    parlo::T
    parhi::T
    perphi::T
end
function CoupledVDF(
    f0; parlower, parupper, perpupper, dpar=nothing, dperp=nothing, normalize=true
)
    plo, phi = promote(float(parlower), float(parupper))
    qhi = oftype(phi, perpupper)
    n = normalize ?
        2π * QuadGK.quadgk(
        v -> v * QuadGK.quadgk(u -> f0(u, v), plo, phi; rtol=1.0e-9)[1],
        zero(qhi), qhi; rtol=1.0e-9
    )[1] : one(plo)
    fn = (u, v) -> f0(u, v) / n
    dp = isnothing(dpar) ? ((u, v) -> _dwrt1(fn, u, v)) : ((u, v) -> dpar(u, v) / n)
    dq = isnothing(dperp) ? ((u, v) -> _dwrt2(fn, u, v)) : ((u, v) -> dperp(u, v) / n)
    return CoupledVDF(fn, dp, dq, plo, phi, qhi)
end

# Regime trait picks the coordinate:
#   NonRelativistic — (p∥,p⊥), pole ζ=(ω−nΩ)/k∥ fixed; outer ∫dp⊥.
#   Relativistic    — (γ,p∥),  pole p∥=(γω−nΩ)/k∥;     outer ∫dγ.
function contribution(d::CoupledVDF, s, ω, k; closure=HarmonicSum())
    return _coupled_contribution(closure, Regime(s), d, s, complex(float(ω)), k)
end

function _coupled_contribution(::HarmonicSum, ::NonRelativistic, d::CoupledVDF, s, ω, k; norm=x -> maximum(abs, x))
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    a = kperp / Ω
    L, U = d.parlo, d.parhi
    p⊥²_mean = 2π * QuadGK.quadgk(
        v -> v^3 * QuadGK.quadgk(u -> d.f0(u, v), L, U; rtol=1.0e-7)[1],
        zero(d.perphi), d.perphi; rtol=1.0e-7
    )[1]
    nmax = nmax_bessel(a^2 * abs(p⊥²_mean) / 2)
    ns = (-nmax):nmax
    χ = first(QuadGK.quadgk(zero(d.perphi), d.perphi; rtol=1.0e-6, norm) do v
        _coupled_perp(v, ns, d, ω, Ω, kz, a, L, U)
    end)
    return SMatrix{3,3,ComplexF64}((s.Pi2 / ω^2) * χ)
end

# Relativistic (γ,p∥) momentum-space path. Momentum distribution f₀ must be
# evaluable at complex p⊥ (the pole pushes p⊥ off-axis).
# Validated vs Maxwell–Jüttner (Swanson) to ~1e-6 and → bi-Maxwellian as μ→∞.
function _coupled_contribution(::HarmonicSum, ::Relativistic, d::CoupledVDF, s, ω, k)
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    a = kperp / Ω
    γmax = sqrt(1 + max(d.parhi^2, d.parlo^2) + d.perphi^2)
    nmax = nmax_bessel(a^2 * d.perphi^2 / 2)
    f = n -> _coupled_harmonic_rel(n, d, ω, Ω, kz, a, γmax)
    χ = converge(f, 1, 1.0e-6; nmax)
    χ = χ .+ _ee33(_bernstein_rel(d, γmax))
    return SMatrix{3,3,ComplexF64}((s.Pi2 / ω^2) * χ)
end

# Relativistic non-resonant e∥e∥ Bernstein term 𝒳_B (derivation §5)
function _bernstein_rel(d, γmax; GLγ=_GLγ, GLp=_GLp)
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
            inner += pw[ip] * (π / 2) * ComplexF64(w * u * d.dpar(u, w) - u^2 * d.dperp(u, w))
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
@inline _U_cov(d, u, w, γ, ω, kz) = kz * d.dpar(u, w) + d.dperp(u, w) * (ω * γ - kz * u) / w

# 3×3 relativistic harmonic integrand 2π·𝒰·𝓣_n at (γ,p∥); bare momenta make 𝓣_n
# regular at w=0. Caller passes w=√(γ²−1−u²) (complex off the real p∥ range).
@inline _rel_integrand(u, w, γ, n, a, ω, kz, d) = (2π * _U_cov(d, u, w, γ, ω, kz)) .* _T_n_bare(n, a * w, u, w)
@inline _rel_integrand(u, γ, n, a, ω, kz, d) = _rel_integrand(u, sqrt(complex(γ^2 - 1 - u^2)), γ, n, a, ω, kz, d)


# One relativistic harmonic, edge-mapped (derivation §5.2.2). 
# Map the disk (γ,p∥) → fixed box (q,θ)∈[0,1]×[−π/2,π/2]:
#   p∥=umax·sinθ, p⊥=umax·cosθ  — inner Jacobian p⊥ cancels the rim 1/p⊥ exactly;
#   γ=1+(γmax−1)q²              — outer Jacobian ∝q flattens the √(γ−1) floor.
# Bessel stays on the fast real path. 
# Off-disk poles (this n doesn't resonate at this γ) aren't peeled — nζ=0 there, so the subtraction reduces to direct integration
function _coupled_harmonic_rel(n, d, ω, Ω, kz, a, γmax; GLγ=_GLγ, GLp=_GLp)
    gn, gw = GLγ
    pn, pw = GLp
    acc = zero(SMatrix{3,3,ComplexF64})
    for ig in eachindex(gn)
        q = (gn[ig] + 1) / 2
        γ = 1 + (γmax - 1) * q^2
        wγ = gw[ig] * (γmax - 1) * q             # gw·½·2(γmax−1)q
        umax = sqrt(γ^2 - 1)
        ζ = (γ * ω - n * Ω) / kz                 # single Landau pole in p∥
        inrange = -umax < real(ζ) < umax
        nζ = inrange ? _rel_integrand(ζ, γ, n, a, ω, kz, d) : zero(SMatrix{3,3,ComplexF64})
        inner = zero(SMatrix{3,3,ComplexF64})
        for ip in eachindex(pn)
            θ = pn[ip] * (π / 2)
            u, w = umax * sin(θ), umax * cos(θ)  # p⊥=w real on the disk
            wu = pw[ip] * (π / 2) * w             # Jacobian p⊥·dθ cancels rim 1/p⊥
            inner = inner .+ wu .* ((_rel_integrand(u, w, γ, n, a, ω, kz, d) .- nζ) ./ (u - ζ))
        end
        inrange && (inner = inner .+ nζ .* _landau_logfac(ζ, -umax, umax))
        acc = acc .+ wγ .* ((-1 / kz) .* inner)
    end
    return acc
end

# I(p⊥) for the WHOLE harmonic sum at one perp node
function _coupled_perp(v, ns, d::CoupledVDF, ω, Ω, kz, a, L, U)
    # Landau–Hilbert for 5 parallel moments: [∂⊥, u·∂⊥, u²·∂⊥, ∂∥, u·∂∥]
    g5(u) = (q=d.dperp(u, v); p=d.dpar(u, v); SVector(q, u * q, u^2 * q, p, u * p))
    ζs = [(ω - n * Ω) / kz for n in ns]
    gζs = g5.(ζs)
    bs = _perp_Bessel_triplet.(ns, a, v)
    # regularized integral part: Σ_n χ_n with the Plemelj removable singularity
    reg = first(QuadGK.quadgk(L, U; rtol=1.0e-7, norm=x -> maximum(abs, x)) do u
        g = g5(u)
        acc = zero(SMatrix{3,3,ComplexF64})
        @inbounds for i in eachindex(ns)
            m = (-1 / kz) .* ((g - gζs[i]) / (u - ζs[i]))
            acc += _In_block((m[1], m[2], m[3], m[4], m[5]), bs[i], v, ω, kz, ns[i] * Ω)
        end
        acc
    end)
    # analytic log-ratio (+ Landau) part, constant in u
    logacc = zero(SMatrix{3,3,ComplexF64})
    @inbounds for i in eachindex(ns)
        m = (-1 / kz) .* (gζs[i] .* _landau_logfac(ζs[i], L, U))
        logacc += _In_block((m[1], m[2], m[3], m[4], m[5]), bs[i], v, ω, kz, ns[i] * Ω)
    end
    return reg + logacc
end

include("qin.jl")