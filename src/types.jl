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
    Species(; Omega, Pi2, vdf, regime=NonRelativistic())

Dimensionless species. `Omega` = signed gyrofrequency ratio
`(Z_s/Z_ref)(m_ref/m_s)`; `Pi2` = `(omega_ps/Omega_ref)^2`.
"""
Base.@kwdef struct Species{T, V, R <: Regime}
    Omega::T          # signed Omega_s_tilde
    Pi2::T            # Pi_s_tilde^2
    vdf::V
    regime::R = Regime(vdf)
end
Species(Omega, Pi2, vdf; regime = Regime(vdf)) = Species(; Omega, Pi2, vdf, regime)

Regime(s::Species) = s.regime

"""
    Plasma(species...)

Container summed over species. Dimensionless throughout; species already hold
`Omega`/`Pi2`, so this is a thin wrapper plus convenience builders.
"""
struct Plasma{S}
    species::S
end
Plasma(species...) = Plasma(Tuple(species))
Plasma(s::Species) = Plasma((s,))

Base.iterate(p::Plasma, state = 1) = iterate(p.species, state)
