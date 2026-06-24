"""
    MaxwellJuttner(; mu)

Relativistic isotropic Maxwell-Juttner VDF. `mu = mc^2 / T` is inverse
dimensionless temperature. Large `mu` tends to `Maxwellian(sqrt(2 / mu))`.
"""
Base.@kwdef struct MaxwellJuttner{T} <: AbstractVDF
    mu::T
end

regime(::MaxwellJuttner) = Relativistic()

# Swanson time-integral form avoids harmonic sums and relativistic resonance
# bookkeeping for isotropic Maxwell-Juttner; ported from LMV.
function contribution(d::MaxwellJuttner, s, ω, k; kwargs...)
    μ = d.mu
    Ω = s.Omega
    kz = para(k)
    kperp = perp(k)
    ω = complex(float(ω))

    if imag(ω) < 0 && iszero(kz)
        throw(ArgumentError("MaxwellJuttner with imag(ω)<0 and kz=0 needs Landau contour continuation"))
    end

    invK2μ = inv(besselkx(2, μ))
    igrand = _maxwell_juttner_swanson_integral(μ, invK2μ, ω, Ω, kz, kperp)
    return SMatrix{3, 3, ComplexF64}((s.Pi2 / ω^2) * (im * ω / Ω * μ^2) * igrand)
end


function _maxwell_juttner_swanson_integral(μ, invK2μ, ω, Ω, kz, kperp)
    σ = sign(real(Ω))
    integrand = ξ -> _maxwell_juttner_integrand(ξ, μ, invK2μ, ω, Ω, kz, kperp)
    return QuadGK.quadgk(
        integrand, 0.0, σ * Inf; rtol = 1.0e-8, atol = 0.0,
        norm = x -> maximum(abs, x)
    )[1]
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
