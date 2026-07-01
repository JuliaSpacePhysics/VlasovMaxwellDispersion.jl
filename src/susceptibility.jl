"""
    contribution(species/vdf, ω, k)

Susceptibility χ_s(ω,k) from one normalized species or vdf.
"""
@inline contribution(s, ω, k; kwargs...) = contribution(s.vdf, s, ω, k; kwargs...)

function contribution(vdf::AbstractVDF, ω, k; kw...)
    return contribution(NormalizedSpecies(1.0, 1.0, vdf), ω, k; kw...)
end


"""
    dielectric(plasma, ω, k; closure=HarmonicSum())

Dielectric tensor `ε = I + Σ_s χ_s(ω,k)`.
"""
function dielectric(plasma, ω, k; kwargs...)
    χ = mapreduce(s -> contribution(s, ω, k; kwargs...), +, NormalizedPlasma(plasma))
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
#   P∂  ← ∫(Bessel)f⊥′    pairs with the f∥ moments M_F^m  (∂⊥ slice)
#   PF  ← ∫(Bessel)f⊥·p⊥  pairs with the f∥′ moments M_T^m  (∂∥ slice)
@inline function _chi_mblock(M, P∂, PF, ω, kz, nΩ)
    MF0, MF1, MF2, MT0, MT1 = M
    # Parallel Landau weights D_m = ω M_F^m − k∥ M_F^{m+1} (∂⊥ slice) and k∥ M_T^m (∂∥ slice).
    # Each tensor entry = (∂⊥ perp bilinear)·wF + (∂∥ perp bilinear)·wT, at order m =
    wF0, wT0 = ω * MF0 - kz * MF1, kz * MT0
    wF1, wT1 = ω * MF1 - kz * MF2, kz * MT1
    xx = P∂[1, 1] * wF0 + PF[1, 1] * wT0
    xy = im * (P∂[1, 2] * wF0 + PF[1, 2] * wT0)
    yy = P∂[2, 2] * wF0 + PF[2, 2] * wT0
    xz = P∂[1, 3] * wF1 + PF[1, 3] * wT1
    zy = im * (P∂[2, 3] * wF1 + PF[2, 3] * wT1)
    zz = nΩ * P∂[3, 3] * MF2 + (ω - nΩ) * PF[3, 3] * MT1   # + non-resonant term
    return SA[xx, xy, xz, yy, zy, zz]
end

# Pointwise (Grid): the perp tensor at node p⊥ before parallel integration
# M=(q,uq,u²q,p,up); here M=c·Δm
@inline function _In_block(Δm, c, bvec, px, ω, kz, nΩ)
    b11, b12, b22, b13, b23, b33 = bvec
    Δ0, Δ1, Δ2, Δ3, Δ4 = Δm
    c2 = 2π * c
    kzpx = kz * px
    D0 = c2 * (ω * Δ0 - kz * Δ1 + kzpx * Δ3)
    D1 = c2 * (ω * Δ1 - kz * Δ2 + kzpx * Δ4)
    zz = (c2 * b33) * (nΩ * Δ2 + (ω - nΩ) * px * Δ4)
    xx, xy, yy = b11 * D0, im * b12 * D0, b22 * D0
    xz, zy = b13 * D1, im * b23 * D1
    return SA[xx, xy, xz, yy, zy, zz]
end

# Materialize the antisymmetric-paire
@inline _antisymmat(t) =
    @SMatrix [t[1] t[2] t[3]; -t[2] t[4] -t[5]; t[3] t[5] t[6]]

# Symmetric 3×3 from its 6 distinct entries (row-major upper triangle).
@inline _symmat(a11, a12, a13, a22, a23, a33) =
    @SMatrix [a11 a12 a13; a12 a22 a23; a13 a23 a33]

# Bessel triplet `bvec = (p⊥Rₙ, p⊥Jₙ′, Jₙ)`
@inline function _perp_Bessel_bilinear(n, a, px)
    z = a * px
    Jm, Jp = besselj(n - 1, z), besselj(n + 1, z)
    b1 = px * (Jm + Jp) / 2
    b2 = px * (Jm - Jp) / 2
    b3 = besselj(n, z)
    return SA[b1 * b1, b1 * b2, b2 * b2, b1 * b3, b2 * b3, b3 * b3]
end

# Whole ±nmax ladder of perp bilinear products
function _perp_Bessel_bilinears(ns, a, px)
    z = a * px
    M = last(ns) + 1
    return @no_escape begin
        Jv = @alloc(typeof(z), M + 1)
        besselj_ladder!(Jv, M, z)
        map(ns) do n
            Jm = _jladder(Jv, n - 1)
            Jp = _jladder(Jv, n + 1)
            Rn = (Jm + Jp) / 2
            Jn = _jladder(Jv, n)
            Jn′ = (Jm - Jp) / 2
            b1 = px * Rn
            b2 = px * Jn′
            SA[b1 * b1, b1 * b2, b2 * b2, b1 * Jn, b2 * Jn, Jn * Jn]
        end
    end
end
