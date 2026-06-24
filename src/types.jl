abstract type AbstractDispersionProblem end

# Orbit-integral evaluator (derivation.md §3), passed as the `closure` keyword of contribution/dielectric/solve.
abstract type IntegralClosure end
struct HarmonicSum <: IntegralClosure end
struct Newberger <: IntegralClosure end

abstract type Regime end
struct NonRelativistic <: Regime end
struct Relativistic <: Regime end

abstract type Continuation end
struct Analytic <: Continuation end       # closed-form Z-function path
struct PiecewisePoly <: Continuation end  # STUB: grid/function VDF

"""
    Wavenumber(kperp, kz)
    Wavenumber(; kz, kperp=zero(kz))

Dimensionless wavevector `k c / Omega_ref`, `(kperp, kz)`,
"""
struct Wavenumber{T}
    kperp::T
    kz::T
end
Wavenumber(kperp, kz) = Wavenumber(promote(kperp, kz)...)
Wavenumber(; kz, kperp = zero(kz)) = Wavenumber(kperp, kz)

@inline para(k::Wavenumber) = k.kz
@inline perp(k::Wavenumber) = k.kperp
@inline Base.abs2(k::Wavenumber) = k.kz^2 + k.kperp^2
@inline Base.angle(k::Wavenumber) = atan(k.kperp, k.kz)  # propagation angle to B0
@inline vec3(k::Wavenumber) = SVector(k.kperp, zero(k.kperp), k.kz)

"""
    NormalizedSpecies(Omega, Pi2, vdf)

Solver's dimensionless per-species representation. `Omega = Ω_s/Ω_ref`; `Pi2 = (ω_ps/Ω_ref)^2`.
"""
Base.@kwdef struct NormalizedSpecies{T, V}
    Omega::T
    Pi2::T
    vdf::V
end

regime(d::NormalizedSpecies) = regime(d.vdf)

"""
    NormalizedPlasma(species...)

Solver's dimensionless container: a bag of [`NormalizedSpecies`](@ref) at one fixed `Ω_ref`.
"""
struct NormalizedPlasma{S} <: AbstractPlasma
    species::S
end
NormalizedPlasma(species::NormalizedSpecies...) = NormalizedPlasma(Tuple(species))
