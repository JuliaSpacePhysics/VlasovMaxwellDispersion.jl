"""
    contribution(species, ω, k::Wavenumber; closure=HarmonicSum())

Dimensionless susceptibility χ_s(ω,k) of one species. Dispatches on the VDF.
"""
@inline contribution(s::Species, ω, k::Wavenumber; kwargs...) = contribution(s.vdf, s, ω, k; kwargs...)


"""
    dielectric(plasma, ω, k; closure=HarmonicSum())

Dielectric tensor `ε = I + Σ_s χ_s(ω,k)`.
"""
function dielectric(plasma, ω, k; kwargs...)
    χ = mapreduce(s -> contribution(s, ω, k; kwargs...), +, plasma)
    return χ + I
end

# Curl-curl operator k̃k̃ᵀ - k̃²I . From the wave eq
# n×(n×E)+εE=0 with n=k̃/ω̃: n×(n×E) = (nnᵀ-n²I)E ⇒ D = ε + curlcurl/ω̃²
@inline function _curlcurl(k)
    kv = vec3(k)
    return kv * kv' - abs2(k) * I
end

"""
    dispersion_tensor(plasma, ω, k::Wavenumber; closure=HarmonicSum())

`𝒟(ω,k) = ε + (k̃k̃ᵀ - k̃²I)/ω̃²`. `det(𝒟)=0` is the dispersion relation.
"""
function dispersion_tensor(plasma, ω, k::Wavenumber; kwargs...)
    ε = dielectric(plasma, ω, k; kwargs...)
    return ε + _curlcurl(k) / complex(float(ω))^2
end

"Aliases for `dispersion_tensor`"
const 𝒟 = dispersion_tensor

"""
    electrostatic_det(plasma, ω, k::Wavenumber) -> ComplexF64

Cheap longitudinal path `k̃ · ε · k̃`; its zeros are the electrostatic modes.
"""
function electrostatic_det(plasma, ω, k::Wavenumber; kwargs...)
    ε = dielectric(plasma, ω, k; kwargs...)
    kv = vec3(k)
    return dot(kv, ε, kv)
end


# --- Shared per-harmonic 3×3 assembler (every VDF) -----------------
# Contracts the cyclotron-harmonic tensor 𝓣_n (derivation §3B) with the parallel
# Landau moments to build χ_n's 3×3 block. The algebra is identical for every
# VDF; only how the moments are obtained differs (Z/Γ_n closed forms for
# Maxwellian vs `hilbert`+Bessel quadrature for arbitrary f). Inputs (§5–6):
#   z = (z0F,z1F,z2F, z0T,z1T)  parallel H∥ moments; F from f∥, T from f∥′
#   p = (JF,J∂F, JdJF,JdJ∂F, ∂J²F,∂J²∂F)  perp Bessel moments (§6 Pⱼ, Pⱼ^∂)
#   nk = nΩ/k⊥  (so nk·k⊥ = nΩ, and z=k⊥p⊥/Ω makes nk a harmonic index over z)
# Every entry is the same Landau combination `D` of one (∂F,F) perp pair with a
# parallel-moment triple; the power of nk indexes the n/z structure of 𝓣_n.
# m33 is the exception: it folds in the non-resonant Bernstein term (§3.2), so it
# uses nΩ and the lower moments M²_F,M¹_T rather than a plain D — see derivation.
@inline function _chi_mblock(z, p, ω, kz, kperp, nk)
    z0F, z1F, z2F, z0T, z1T = z
    D(X, Y, a, b, c) = kz * (-X * a + Y * b) + ω * X * c
    m11 = nk^2 * D(p.J∂F, p.JF, z1F, z0T, z0F)
    m21 = -im * nk * D(p.JdJ∂F, p.JdJF, z1F, z0T, z0F)
    m31 = nk * D(p.J∂F, p.JF, z2F, z1T, z1F)
    m22 = D(p.∂J²∂F, p.∂J²F, z1F, z0T, z0F)
    m32 = im * D(p.JdJ∂F, p.JdJF, z2F, z1T, z1F)
    m33 = kperp * nk * p.J∂F * z2F + (ω - kperp * nk) * p.JF * z1T
    return @SMatrix ComplexF64[m11 -m21 m31; m21 m22 -m32; m31 m32 m33]
end


# The 6 perp Bessel moments at one v⊥ node:
#  2π·{powers of v⊥}·{Jₙ², JₙJₙ′, Jₙ′²} at z=a·v⊥. The f⊥ weight and
# the ∫dp⊥ are supplied by the caller's quadrature.
@inline function _perp_bessel_moments(n, a, v)
    Jn, Jn′ = besselj(n, a * v), _besselj_prime(n, a * v)
    tp = 2π
    return (
        JF=tp * v * Jn^2, J∂F=tp * Jn^2,
        JdJF=tp * v^2 * Jn * Jn′, JdJ∂F=tp * v * Jn * Jn′,
        ∂J²F=tp * v^3 * Jn′^2, ∂J²∂F=tp * v^2 * Jn′^2,
    )
end
