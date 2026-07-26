# Literal shifted-Gaussian ring
# no finite Bessel closure (docs/Maxwellian.md); perp moments use Route A — the exact-shift
# parabolic-cylinder series. Accurate for Λr=k⊥ vr/Ω ≲ 10; beyond, use `SeparableVDF`.

"""
    GaussianRing(; vth, vd=0, vr=0)

Drifting ring-beam with a *literal* shifted-Gaussian perpendicular factor:

    f ∝ Gaussian(vth_perp, vr) ⊗ Gaussian(vth_para, vd)

`vth` is a scalar (isotropic) or a `(⊥, ∥)` pair.
Accurate for `Λr=k⊥ vr/Ω ≲ 10`; beyond, use [`SeparableVDF`](@ref).
`vr=nothing` reduces to the bi-Maxwellian.
"""
function GaussianRing(; vth, vd = nothing, vr = nothing)
    vth_perp, vth_para = perp_para(vth)
    perp = (isnothing(vr) || iszero(vr)) ? Gaussian(vth_perp) : Gaussian(vth_perp, vr)
    return perp ⊗ Gaussian(vth_para, @something(vd, zero(vth_para)))
end


struct GaussianRingCtx{T}
    vr::T
    p::T             # 1/vth²
    P::Vector{T}     # parabolic-cylinder moment table 𝓔ₖ
    Lcap::Int
    nmax::Int
end

# Drifted perp Gaussian → Route-A table (the magnitude ring). `vd` is the perp shift vr.
function perp_setup(d::Gaussian{<:Any, <:Real}, β)
    iszero(d.vd) && return perp_setup(Gaussian(d.vth), β)
    vr, p = promote(d.vd, 1 / d.vth^2)
    mwin = nmax_bessel((β * vr)^2 / 2)                       # ~Λr harmonic/series reach
    nλ = nmax_bessel((d.vth^2 / 2) * β^2)
    nmax = nλ + mwin + 2
    Lcap = nλ + mwin + 8
    P = _paracyl_moments(vr, p, 2 * (nmax + 1) + 2 * Lcap + 4)
    return GaussianRingCtx(vr, p, P, Lcap, nmax)
end
nmax_harm(c::GaussianRingCtx, β) = c.nmax

# Parabolic-cylinder moments 𝓔ₖ = ∫₀^∞ vᵏ e^{-(v-vr)²p} dv (P[k+1]=𝓔ₖ, p=1/vth²),
# erfc-seeded two-term recurrence 𝓔ₖ = vr·𝓔_{k-1} + (k-1)/(2p)·𝓔_{k-2}.
function _paracyl_moments(vr, p, kmax)
    P = Vector{typeof(float(vr))}(undef, kmax + 1)
    P[1] = sqrt(π / p) * erfc(-vr * sqrt(p)) / 2
    P[2] = vr * P[1] + exp(-p * vr^2) / (2p)
    for k in 2:kmax
        P[k + 1] = vr * P[k] + (k - 1) / (2p) * P[k - 1]
    end
    return P
end

# Route-A moment  𝒮_q = ∫₀^∞ v^q e^{-(v-vr)²p} J_μ(βv) J_ν(βv) dv  via the Bessel-product
# power series  J_μJ_ν=Σ_l e^{μν}_l (βv)^{μ+ν+2l}, integrated term-by-term against P.
# Negative orders fold in by J_{-m}=(-1)^m J_m. `e^{μν}_0=1/(μ!ν!2^{μ+ν})`.
@inline function _jj_moment(q, μ, ν, β, P, Lcap, tol)
    sgn = 1
    μ < 0 && (isodd(μ) && (sgn = -sgn); μ = -μ)
    ν < 0 && (isodd(ν) && (sgn = -sgn); ν = -ν)
    base = μ + ν
    q + base + 1 > length(P) && return zero(eltype(P))
    term = exp(-loggamma(μ + 1) - loggamma(ν + 1) - base * log(2.0)) * β^base * P[q + base + 1]
    s = term
    for l in 0:(Lcap - 1)
        q + base + 2l + 3 > length(P) && break
        r = -(base + 2l + 2) * (base + 2l + 1) /
            (4 * (l + 1) * (base + l + 1) * (μ + l + 1) * (ν + l + 1))
        term *= r * β^2 * P[q + base + 2l + 3] / P[q + base + 2l + 1]
        s += term
        abs(term) <= tol * abs(s) && break
    end
    return sgn * s
end

# Perp tensor from the gyro polarization triplet (p⊥Rₙ, p⊥Jₙ′, Jₙ), with p⊥Rₙ kept as the
# moment ⟨p⊥(J_{n-1}+J_{n+1})/2⟩ instead of (n/β)⟨Jₙ⟩ — so every entry is a genuine velocity
# moment of {J_{n-1},J_n,J_{n+1}} products, finite at β=0 (no n/β). The Rₙ-rows reuse the
# same moments as the Jₙ′ entries, with the cross-order sign flipped (+ for Rₙ², − for Jₙ′²).
# `g(q,μ,ν)` is the q-th p⊥-moment of J_μJ_ν; `q0` is the slice's base p⊥ power.
@inline function _ring_perp_tensor(g, n, q0)
    Jmm = g(q0 + 2, n - 1, n - 1)
    Jpp = g(q0 + 2, n + 1, n + 1)
    Jmp = g(q0 + 2, n - 1, n + 1)
    Sm = g(q0 + 1, n, n - 1)
    Sp = g(q0 + 1, n, n + 1)
    J² = g(q0, n, n)                       # Jₙ² · [3,3]
    JdJ = (Sm - Sp) / 2                    # p⊥JₙJₙ′ · [2,3]
    Jd² = (Jmm - 2Jmp + Jpp) / 4           # p⊥²Jₙ′² · [2,2]
    RJ = (Sm + Sp) / 2                     # p⊥Rₙ·Jₙ · [1,3]
    RJd = (Jmm - Jpp) / 4                  # p⊥Rₙ·p⊥Jₙ′ · [1,2]
    R² = (Jmm + 2Jmp + Jpp) / 4            # (p⊥Rₙ)² · [1,1]
    return _symmat(R², RJd, RJ, Jd², JdJ, J²)
end

@inline function perp_moments(c::GaussianRingCtx, n, β)
    vr, p, P, Lcap = c.vr, c.p, c.P, c.Lcap
    tol = oftype(p, 1.0e-8)
    sm(q, μ, ν) = _jj_moment(q, μ, ν, β, P, Lcap, tol) / P[2]
    sd(q, μ, ν) = 2p * (vr * sm(q, μ, ν) - sm(q + 1, μ, ν))    # ∂F slice via f⊥′
    P∂ = _ring_perp_tensor(sd, n, 0)                          # ∂F slice: base power q0=0
    PF = _ring_perp_tensor(sm, n, 1)                          # F slice : base power q0=1
    return P∂, PF
end
