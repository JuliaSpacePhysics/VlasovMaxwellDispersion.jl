"""
    Gaussian(vth, vd=nothing)

Drifting Gaussian 1D factor `∝ exp(-(v-vd)²/vth²)`.
"""
struct Gaussian{T, D}
    vth::T
    vd::D
end

Gaussian(vth) = Gaussian(vth, nothing)
parallel_even(g::Gaussian) = g.vd === nothing || iszero(g.vd)

# unnormalized 1d density shape
@inline (g::Gaussian)(v) = exp(-((v - @something(g.vd, zero(g.vth))) / g.vth)^2)

# `Z(z) = i√π · w(z)` where Faddeeva function `w(z) = erfcx(-i z)`.
@inline plasma_dispersion_function(z) = im * sqrt(oftype(real(z), pi)) * erfcx(-im * z)

const Z = plasma_dispersion_function

function para_moments(p::Gaussian, Δ, kz)
    vthpar = p.vth
    vd = @something p.vd zero(p.vth)
    if iszero(kz)
        # M_F^m = ⟨uᵐ⟩/Δ, M_T^m = ∫uᵐf′/Δ
        invΔ = 1 / Δ
        return (invΔ, vd * invΔ, (vd^2 + vthpar^2 / 2) * invΔ, zero(invΔ), -invΔ)
    end
    σ⁻¹ = 1 / (kz * vthpar)
    ζ = (Δ - kz * vd) * σ⁻¹
    # k∥<0 flips the Landau contour above the pole
    Z0 = kz > 0 ? Z(ζ) : -Z(-ζ)
    Z1 = 1 + ζ * Z0
    Z2 = ζ * Z1
    MF0 = -Z0 * σ⁻¹
    MF1 = -(Z0 * vd + Z1 * vthpar) * σ⁻¹
    MF2 = -(Z0 * vd^2 + Z1 * 2 * vthpar * vd + Z2 * vthpar^2) * σ⁻¹
    invth2 = 2 / vthpar^2
    MT0 = (MF0 * vd - MF1) * invth2
    MT1 = (MF1 * vd - MF2) * invth2
    return (MF0, MF1, MF2, MT0, MT1)
end
