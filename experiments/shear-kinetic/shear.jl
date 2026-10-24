# Local uniform-shear electrostatic dielectric: theory.md §3 eq. (10). Scope: §5.
# Package units — ω and S normalized to Ω_ref, speeds in c, k in Ω_ref/c. A VDF handed to
# NormalizedSpecies is read in the invariant w of theory.md eq. (5).

using VlasovMaxwellDispersion
using LinearAlgebra: dot
using StaticArrays: SVector

"""Wavevector with an explicit `kₓ`/`kᵧ` split — shear breaks the perpendicular rotational
symmetry that lets `Wavenumber` carry only `(k⊥,k∥)`."""
struct ShearK{T}
    kx::T
    ky::T
    kz::T
end
ShearK(; kx = 0.0, ky, kz) = ShearK(promote(kx, ky, kz)...)
Base.abs2(k::ShearK) = k.kx^2 + k.ky^2 + k.kz^2

"Species in the circularizing coordinates of theory.md §1; returns `(species′, k_eff)`."
function shear_map(s::NormalizedSpecies, k::ShearK, S)
    f = 1 + S / s.Omega                       # 1 + σ_s
    f > 0 || throw(DomainError(f, "σ ≤ −1: orbits unbound, local theory invalid"))
    keff = sqrt(k.kx^2 + k.ky^2 / f)
    return NormalizedSpecies(s.Omega * sqrt(f), s.Pi2, s.vdf), keff
end

"""Longitudinal `ε_L(ω,k)` under linear shear `S`. `u0` is the E×B speed at the guiding
centre the analysis is local to — a frame choice, kept explicit. Neglects theory.md eq. (12);
error ≈ 1.2σ."""
function shear_epsilon(plasma, ω, k::ShearK, S; u0 = 0.0, kw...)
    k2 = abs2(k)
    ε = one(complex(float(real(ω))))
    for s in NormalizedPlasma(plasma).species
        s′, keff = shear_map(s, k, S)
        kv = SVector(keff, zero(keff), k.kz)
        χ = contribution(s′, ω - k.ky * u0, Wavenumber(keff, k.kz); kw...)
        ε += dot(kv, χ, kv) / k2
    end
    return ε
end

"Unsheared `k̂·ε·k̂` at the same physical `k`; equals `shear_epsilon(…, S=0)`."
function epsilon_L(plasma, ω, k::ShearK; kw...)
    kperp = sqrt(k.kx^2 + k.ky^2)
    kv = SVector(kperp, zero(kperp), k.kz)
    ε = one(complex(float(real(ω))))
    for s in NormalizedPlasma(plasma).species
        ε += dot(kv, contribution(s, ω, Wavenumber(kperp, k.kz); kw...), kv) / abs2(k)
    end
    return ε
end
