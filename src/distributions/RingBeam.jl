# Literal shifted-Gaussian ring-beam: perp factor is a Gaussian in the *magnitude*
# v⊥ (eq.(13) form), NOT the gyro-averaged I₀ ring of `Maxwellian(; vr)`. This form
# has no finite Bessel closure (docs/Maxwellian.md); perp moments use Route A —
# the exact-shift parabolic-cylinder series. Parallel factor is the drifting
# Gaussian, so its Z-moments are shared verbatim with the bi-Maxwellian.

"""
    RingBeam(; vth_par, vth_perp=vth_par, vd=0, vr=0)

Drifting **ring-beam** with a *literal* shifted-Gaussian perpendicular factor:

    f ∝ exp(-(v∥-vd)²/vth_par²) · exp(-(v⊥-vr)²/vth_perp²)

The perp factor is the magnitude form `exp(-(v⊥-vr)²/c⊥²)` (normaliser
`A=e^{-(vr/c⊥)²}+√π(vr/c⊥)erfc(-vr/c⊥)`) used in ring-beam instability studies —
distinct from [`Maxwellian`](@ref)`(; vr)`, whose ring is the gyro-averaged `I₀`
form. It has **no finite Bessel closure**; perp moments use the exact-shift
parabolic-cylinder series ("Route A", `docs/Maxwellian.md`): `Jₙ(βv)²` expanded in
powers, integrated against the erfc-seeded moments `𝓔ₖ=∫₀^∞ vᵏ e^{-(v-vr)²/c⊥²}dv`.

Accurate for `Λr=k⊥ vr/Ω ≲ 10`; beyond, the alternating series loses precision
(use [`SeparableVDF`](@ref)). `vr=0` reduces to the bi-Maxwellian fast path.
"""
Base.@kwdef struct RingBeam{T} <: AbstractVDF
    vth_par::T
    vth_perp::T = vth_par
    vd::T = zero(vth_par)
    vr::T = zero(vth_par)
end

@inline thermal_par(d::RingBeam) = d.vth_par
@inline thermal_perp(d::RingBeam) = d.vth_perp
@inline drift(d::RingBeam) = d.vd
@inline ring(d::RingBeam) = d.vr

# Parabolic-cylinder moments 𝓔ₖ = ∫₀^∞ vᵏ e^{-(v-pr)²p} dv (E[k+1]=𝓔ₖ, p=1/p_thperp²),
# erfc-seeded two-term recurrence 𝓔ₖ = pr·𝓔_{k-1} + (k-1)/(2p)·𝓔_{k-2}.
function _paracyl_moments(pr, p, kmax)
    P = Vector{Float64}(undef, kmax + 1)
    P[1] = sqrt(π / p) * erfc(-pr * sqrt(p)) / 2
    P[2] = pr * P[1] + exp(-p * pr^2) / (2p)
    for k in 2:kmax
        P[k+1] = pr * P[k] + (k - 1) / (2p) * P[k-1]
    end
    return P
end

# Route-A moment  𝒮_q = ∫₀^∞ v^q e^{-(v-pr)²p} J_μ(βv) J_ν(βv) dv  via the Bessel-product
# power series  J_μJ_ν=Σ_l e^{μν}_l (βv)^{μ+ν+2l}, integrated term-by-term against P.
# Negative orders fold in by J_{-m}=(-1)^m J_m. `e^{μν}_0=1/(μ!ν!2^{μ+ν})`.
@inline function _jj_moment(q, μ, ν, β, P, Lcap, tol)
    sgn = 1
    μ < 0 && (isodd(μ) && (sgn = -sgn); μ = -μ)
    ν < 0 && (isodd(ν) && (sgn = -sgn); ν = -ν)
    base = μ + ν
    q + base + 1 > length(P) && return zero(eltype(P))
    term = exp(-loggamma(μ + 1) - loggamma(ν + 1) - base * log(2.0)) * β^base * P[q+base+1]
    s = term
    for l in 0:Lcap-1
        q + base + 2l + 3 > length(P) && break
        r = -(base + 2l + 2) * (base + 2l + 1) /
            (4 * (l + 1) * (base + l + 1) * (μ + l + 1) * (ν + l + 1))
        term *= r * β^2 * P[q+base+2l+3] / P[q+base+2l+1]
        s += term
        abs(term) <= tol * abs(s) && break
    end
    return sgn * s
end

function contribution(d::RingBeam, s, ω, k; rtol=1.0e-8, kwargs...)
    Ω = s.Omega
    kz = para(k)
    kperp = perp(k)
    cperp = thermal_perp(d)
    ω = complex(float(ω))
    β = kperp / Ω
    λ = (cperp^2 / 2) * β^2
    prefac = s.Pi2 / ω^2
    if iszero(d.vr) || iszero(kperp)                       # plain bi-Maxwellian fast path
        nmax = nmax_bessel(λ)
        f = n -> _maxwellian_harmonic(n, ω, Ω, kz, β, d.vth_par, cperp, d.vd)
        return prefac * converge(f, 1, rtol; nmax)
    end
    pr = d.vr                                              # perp ring momentum
    p = 1 / cperp^2
    mwin = nmax_bessel((β * pr)^2 / 2)                     # ~Λr harmonic/series reach
    nmax = nmax_bessel(λ) + mwin + 2
    Lcap = nmax_bessel(λ) + mwin + 8
    P = _paracyl_moments(pr, p, 2 * (nmax + 1) + 2 * Lcap + 4)
    Pnorm = P[2]                                           # 𝓔₁ = ∫v e^{-(v-pr)²p}dv (perp norm)
    f = n -> _ringbeam_harmonic(n, ω, Ω, kz, β, d.vth_par, d.vd, pr, p, P, Pnorm, Lcap, rtol)
    return prefac * converge(f, 1, rtol; nmax)
end

# One cyclotron harmonic. Parallel Z-moments identical to the bi-Maxwellian; perp
# moments are the literal-ring Route-A fundamentals (Jₙ², v⊥JₙJₙ′, v⊥²Jₙ′²) in the
# F-slice (×v⊥f⊥) and ∂F-slice (×f⊥′, with f⊥′=-2p(v-pr)e^{-(v-pr)²p}).
@inline function _ringbeam_harmonic(n, ω, Ω, kz, β, cpar, vd, pr, p, P, Pnorm, Lcap, tol)
    nΩ = n * Ω
    M = _gaussian_par_moments(ω, kz, nΩ, cpar, vd)

    sm(q, μ, ν) = _jj_moment(q, μ, ν, β, P, Lcap, tol) / Pnorm
    sd(q, μ, ν) = 2p * (pr * sm(q, μ, ν) - sm(q + 1, μ, ν))    # ∂F slice via f⊥′
    AF = sm(1, n, n)
    BF = (sm(2, n, n - 1) - sm(2, n, n + 1)) / 2
    CF = (sm(3, n - 1, n - 1) - 2 * sm(3, n - 1, n + 1) + sm(3, n + 1, n + 1)) / 4
    A∂ = sd(0, n, n)
    B∂ = (sd(1, n, n - 1) - sd(1, n, n + 1)) / 2
    C∂ = (sd(2, n - 1, n - 1) - 2 * sd(2, n - 1, n + 1) + sd(2, n + 1, n + 1)) / 4

    nβ = n / β
    P∂ = _symmat(nβ^2 * A∂, nβ * B∂, nβ * A∂, C∂, B∂, A∂)
    PF = _symmat(nβ^2 * AF, nβ * BF, nβ * AF, CF, BF, AF)
    return _chi_mblock(M, P∂, PF, ω, kz, nΩ)
end
