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
    χ = mapreduce(s -> contribution(s, ω, k; kwargs...), +, Plasma(plasma))
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


# Builds one cyclotron-harmonic block χ_n by contracting the perp Bessel tensor
# with the parallel Landau moments (derivation §5.1). Same algebra for every VDF;
# only how the moments are obtained differs (Z/Γ_n closed forms for Maxwellian vs
# `hilbert`+Bessel quadrature for arbitrary f).
#
# The numerator p⊥U splits into a ∂f/∂p⊥ and a ∂f/∂p∥ gradient slice, giving two
# perp Bessel-bilinear matrices and two parallel-moment families:
#   P∂  ← ∫(Bessel)f⊥′    pairs with the f∥ moments z*F   (the M_F / ∂⊥ slice)
#   PF  ← ∫(Bessel)f⊥·p⊥  pairs with the f∥′ moments z*T   (the M_T / ∂∥ slice)
@inline function _chi_mblock(z, P∂, PF, ω, kz, nΩ)
    z0F, z1F, z2F, z0T, z1T = z
    # Parallel Landau weights D_m = ω M_F^m − k∥ M_F^{m+1} (∂⊥ slice) and k∥ M_T^m (∂∥ slice).
    # Each tensor entry = (∂⊥ perp bilinear)·wF + (∂∥ perp bilinear)·wT, at order m =
    wF0, wT0 = ω * z0F - kz * z1F, kz * z0T
    wF1, wT1 = ω * z1F - kz * z2F, kz * z1T
    xx = P∂[1, 1] * wF0 + PF[1, 1] * wT0
    xy = im * (P∂[1, 2] * wF0 + PF[1, 2] * wT0)
    yy = P∂[2, 2] * wF0 + PF[2, 2] * wT0
    xz = P∂[1, 3] * wF1 + PF[1, 3] * wT1
    yz = im * (P∂[2, 3] * wF1 + PF[2, 3] * wT1)
    zz = nΩ * P∂[3, 3] * z2F + (ω - nΩ) * PF[3, 3] * z1T   # + non-resonant term
    return @SMatrix ComplexF64[xx xy xz; -xy yy -yz; xz yz zz]
end

# Pointwise (Coupled/Grid): the perp tensor at node v⊥ before parallel integration
@inline function _In_block(z, bvec, px, ω, kz, nΩ)
    b1, b2, b3 = bvec
    z0F, z1F, z2F, z0T, z1T = z
    D0 = 2π * (ω * z0F - kz * z1F + kz * px * z0T)
    D1 = 2π * (ω * z1F - kz * z2F + kz * px * z1T)
    xx, xy, yy = b1 * b1 * D0, im * b1 * b2 * D0, b2 * b2 * D0
    xz, yz = b1 * b3 * D1, im * b2 * b3 * D1
    zz = 2π * b3 * b3 * (nΩ * z2F + (ω - nΩ) * px * z1T)
    return @SMatrix ComplexF64[xx xy xz; -xy yy -yz; xz yz zz]
end

# Symmetric 3×3 from its 6 distinct entries (row-major upper triangle).
@inline _symmat(a11, a12, a13, a22, a23, a33) =
    @SMatrix [a11 a12 a13; a12 a22 a23; a13 a23 a33]

# Bessel triplet `bvec = (p⊥Rₙ, p⊥Jₙ′, Jₙ)`
@inline function _perp_Bessel_triplet(n, a, px)
    z = a * px
    Jm, Jp = besselj(n - 1, z), besselj(n + 1, z)
    return SVector(px * (Jm + Jp) / 2, px * (Jm - Jp) / 2, besselj(n, z))
end
