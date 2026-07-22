"""
    MaxwellJuttner(; mu)

Relativistic isotropic Maxwell-Juttner VDF. `mu = mc^2 / T` is inverse
dimensionless temperature. Large `mu` tends to `Maxwellian(sqrt(2 / mu))`.
"""
Base.@kwdef struct MaxwellJuttner{T} <: AbstractVDF
    mu::T
end

regime(::MaxwellJuttner) = Relativistic()
parallel_even(::MaxwellJuttner) = true

# Relativistic density f(q,u)∝exp(-μγ), γ=√(1+q²+u²). Feeds the general CoupledVDF path.
@inline (d::MaxwellJuttner)(q, u) = exp(-d.mu * sqrt(1 + q^2 + u^2))

# Swanson time-integral form avoids harmonic sums and relativistic resonance
# bookkeeping for isotropic Maxwell-Juttner; ported from LMV.
function contribution(d::MaxwellJuttner, s, ω, k; kwargs...)
    μ = d.mu
    Ω = s.Omega
    kz = para(k)
    kperp = perp(k)

    if imag(ω) < 0 && iszero(kz)
        throw(ArgumentError("MaxwellJuttner with imag(ω)<0 and kz=0 needs Landau contour continuation"))
    end
    invK2μ = inv(besselkx(2, μ))
    # near-marginal layer (|Im ω| ≤ 1e-4|Ω|) stays on the straight integral: it
    # approximates the physical boundary value to O(Im ω), where the continued
    # sheet is exponentially far away (see _mj_parallel_superluminal_integral)
    igrand = if imag(ω) < -1.0e-4 * abs(Ω) && real(ω)^2 > kz^2
        if iszero(kperp)
            _mj_parallel_superluminal_integral(μ, invK2μ, ω, Ω, kz)
        else
            return _mj_cycle_contribution(d, s, ω, k; kwargs...)
        end
    else
        _maxwell_juttner_swanson_integral(μ, invK2μ, ω, Ω, kz, kperp)
    end
    return SMatrix{3, 3, ComplexF64}((s.Pi2 / ω^2) * (im * ω / Ω * μ^2) * igrand)
end

# maxevals bounds the damped-superluminal case (integrand stops decaying,
# quadgk would grind unboundedly on solver iterates passing through)
function _maxwell_juttner_swanson_integral(μ, invK2μ, ω, Ω, kz, kperp; rtol = 1.0e-8, norm = NORM, maxevals = 100_000)
    σ = sign(real(Ω))
    integrand = ξ -> _maxwell_juttner_integrand(ξ, μ, invK2μ, ω, Ω, kz, kperp)
    return QuadGK.quadgk(integrand, 0.0, σ * Inf; rtol, norm, maxevals)[1]
end

function _maxwell_juttner_integrand(ξ, μ, invK2μ, ω, Ω, kz, kperp)
    sinξ, cosξ = sincos(ξ)
    T1 = @SMatrix ComplexF64[
        cosξ sinξ 0;
        -sinξ cosξ 0;
        0 0 1
    ]
    Qxx = kperp^2 * sinξ^2
    Qxy = kperp^2 * sinξ * (1 - cosξ)
    Qxz = kperp * kz * ξ * sinξ
    Qyy = -kperp^2 * (1 - cosξ)^2
    Qyz = -kperp * kz * ξ * (1 - cosξ)
    Qzz = kz^2 * ξ^2
    T2 = (1 / Ω)^2 * @SMatrix ComplexF64[
        Qxx Qxy Qxz;
        -Qxy Qyy Qyz;
        Qxz -Qyz Qzz
    ]
    R = ((μ * Ω - im * ξ * ω)^2 + 2 * kperp^2 * (1 - cosξ) + (kz * ξ)^2) / Ω^2
    sqrtR = sqrt(R)
    real(sqrtR) < 0 && (sqrtR = -sqrtR)
    return (
        _besselk_ratio(2, sqrtR, μ, invK2μ) * T1 -
            _besselk_ratio(3, sqrtR, μ, invK2μ) / sqrtR * T2
    ) / R
end


@inline function _besselk_ratio(ν, sqrtR, μ, invK2μ)
    Kν = abs(sqrtR) > 1.0e6 ? sqrt(π / (2sqrtR)) : besselkx(ν, sqrtR)
    return exp(μ - sqrtR) * Kν * invK2μ
end

# At parallel damped-superluminal k, the straight ξ-integral is not the
# analytic continuation (the branch continued from the subluminal side diverges
# on the real path). With ξ = σs, R factorizes as R₀ = A·B/Ω²,
#   A = μ|Ω| − is(ω−k∥),  B = μ|Ω| − is(ω+k∥),
# so √R₀ = √A·√B with per-factor principal sqrts is the continued branch along
# any path avoiding the two branch points s*∓ = −iμ|Ω|/(ω∓k∥), whose cuts run
# radially outward. The subluminal-germ continuation integrates along a ray in
# the wedge between the two cut rays (the same side of each branch point as the
# real path had before ω crossed the light line). Harmonic components e^{ijs}
# are split — folded into exp(μ−√R₀+ijs) so growth and decay combine
# analytically — because the cyclotron-resonant ones can lack a convergent
# in-wedge ray (shallow damping): those take a dogleg that bends to a
# convergent direction after clearing the branch-point radii, flipping the
# √ branch at each cut crossing. Oblique continuation instead uses transported
# residue cycles. Near-marginal in-band ω (Im ω → 0⁻) is
# exponentially far from the UHP boundary value on this sheet — evaluate such
# roots at Im ω ≥ 0 instead.
# ---- Transported residue cycles (every k⊥ ≠ 0 damped-superluminal point routes here) ----
# Thin specialization of the generic machinery (residue_cycles.jl): the MJ
# energy-form gradient is (∂_γf, ∂_∥f) = (−μe^{−μγ}, 0), so 𝒰 = −μω e^{−μγ}.

const _MJ_BOX_CACHE = Dict{Float64, Any}()
const _MJ_CYCLE_QUAD = BoxQuad(GaussLegendre(40), GaussLegendre(48))

# (q,u) box for the straight part: e^{−μ(γmax−1)} ≈ 1e-7 truncation (wider boxes
# lose more to fixed-quadrature resolution than they gain in truncation)
_mj_prepared_box(μ) = get!(_MJ_BOX_CACHE, μ) do
    P = sqrt((1 + 16 / μ)^2 - 1)
    prepare(CoupledVDF((γ, u) -> exp(-μ * γ); para = (-P, P), perp = P,
        coords = :energy, regime = Relativistic()); quad = _MJ_CYCLE_QUAD)
end

function _mj_cycle_contribution(d::MaxwellJuttner, s, ω, k; quad = _MJ_CYCLE_QUAD, kwargs...)
    c = _mj_prepared_box(d.mu)
    logU0 = log(complex(-d.mu * ω))
    scaledUcov = (γ, u, σ) -> exp(logU0 - d.mu * γ + σ)
    return _cycle_contribution(c, s, ω, k; mueff = d.mu, scaledUcov, quad, kwargs...) / c.cache.n
end

function _mj_parallel_superluminal_integral(μ, invK2μ, ω, Ω, kz; kw...)
    σ = sign(real(Ω))
    aΩ = abs(Ω)
    Jp = _mj_J(μ, invK2μ, ω, aΩ, kz, 1, 0, 2; rtol = 2.0e-9, kw...)
    Jm = _mj_J(μ, invK2μ, ω, aΩ, kz, -1, 0, 2; rtol = 2.0e-9, kw...)
    J0 = _mj_J(μ, invK2μ, ω, aΩ, kz, 0, 0, 2; rtol = 1.0e-9, kw...)
    bzz = (kz / aΩ)^2
    Jz = _mj_J(μ, invK2μ, ω, aΩ, kz, 0, 2, 3;
        rtol = clamp(1.0e-9 / bzz, 1.0e-9, 1.0e-4), kw...)
    xx = (Jp + Jm) / 2
    xy = σ * im * (Jm - Jp) / 2
    zz = J0 - bzz * Jz
    return σ * @SMatrix ComplexF64[xx xy 0; -xy xx 0; 0 0 zz]
end

function _mj_J_component(s, μ, invK2μ, ω, aΩ, kz, flip, j, q, ν)
    A = μ * aΩ - im * s * (ω - kz)
    B = μ * aΩ - im * s * (ω + kz)
    sqrtR = flip * sqrt(A) * sqrt(B) / aΩ
    R = A * B / aΩ^2
    E = exp(μ - sqrtR + im * j * s) * s^q
    Kx = abs(sqrtR) > 1.0e6 ? sqrt(π / (2sqrtR)) : besselkx(ν, sqrtR)
    return E * Kx * sqrtR^(2 - ν) / R * invK2μ
end

# asymptotic √R direction along the ray s = t·cis(φ): √R ≈ t·ĉ(φ)
_mj_chat(φ, ω, aΩ, kz) =
    sqrt(-im * cis(φ) * (ω - kz)) * sqrt(-im * cis(φ) * (ω + kz)) / aΩ
_mj_rate(φ, ω, aΩ, kz, harm) = real(_mj_chat(φ, ω, aΩ, kz)) + harm * sin(φ)

function _mj_J(μ, invK2μ, ω, aΩ, kz, j, q, ν; rtol = 1.0e-9, maxevals = 10^5)
    ψm = -π / 2 - angle(ω - kz)
    ψp = -π / 2 - angle(ω + kz)
    ψlo, ψhi = minmax(ψm, ψp)
    pad = 0.02 * (ψhi - ψlo)
    cand = range(ψlo + pad, ψhi - pad, length = 7)
    rates = map(φ -> _mj_rate(φ, ω, aΩ, kz, j), cand)
    ib = argmax(rates)
    φ1 = cand[ib]
    e1 = cis(φ1)
    if rates[ib] > 0.05
        return quadgk(t -> _mj_J_component(t * e1, μ, invK2μ, ω, aΩ, kz, 1, j, q, ν) * e1,
            0.0, Inf; rtol, maxevals)[1]
    end
    rm = μ * aΩ / abs(ω - kz)
    rp = μ * aΩ / abs(ω + kz)
    R0 = 2 * max(rm, rp) + 5
    M = R0 * e1
    cand2 = range(ψhi + 0.05, π - 0.05, length = 40)
    rates2 = map(φ -> -real(_mj_chat(φ, ω, aΩ, kz)) + j * sin(φ), cand2)
    i2 = argmax(rates2)
    rates2[i2] > 0.02 || return complex(NaN, NaN)   # no convergent path; NaN-guarded upstream
    e2 = cis(cand2[i2])
    crossings = Float64[]
    for (ψ, rb) in ((ψm, rm), (ψp, rp))
        den = imag(e2 * cis(-ψ))
        abs(den) < 1.0e-12 && continue
        t = -imag(M * cis(-ψ)) / den
        t > 0 || continue
        real((M + t * e2) * cis(-ψ)) ≥ 0.99 * rb && push!(crossings, t)
    end
    sort!(crossings)
    total = quadgk(t -> _mj_J_component(t * e1, μ, invK2μ, ω, aΩ, kz, 1, j, q, ν) * e1,
        0.0, R0; rtol, maxevals)[1]
    knots = vcat(0.0, crossings)
    for (i, t0) in enumerate(knots)
        t1 = i < length(knots) ? knots[i + 1] : Inf
        flip = iseven(i - 1) ? 1 : -1
        total += quadgk(t -> _mj_J_component(M + t * e2, μ, invK2μ, ω, aΩ, kz, flip, j, q, ν) * e2,
            t0, t1; rtol, maxevals)[1]
    end
    return total
end
