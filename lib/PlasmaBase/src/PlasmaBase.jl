"""
    PlasmaBase

Thin, solver-agnostic vocabulary for describing a plasma *physically*.
It holds the PROBLEM description (identity, populations, the system and its field) plus an accessor interface.
"""
module PlasmaBase

export AbstractVDF, Particle, Electron, Proton, Species, Plasma
export charge, mass, particle, number_density, distribution, species, magnetic_field
export gyrofrequency_ratio, plasma_gyro_ratio

abstract type AbstractVDF end
abstract type AbstractPlasma end

# --- SI constants (plain Float64) ---
const C_SI = 2.99792458e8       # m/s
const E_SI = 1.602176634e-19    # C
const ME_SI = 9.1093837015e-31   # kg
const MP_SI = 1.67262192369e-27  # kg
const EPS0_SI = 8.8541878128e-12   # F/m
const KB_SI = 1.380649e-23       # J/K

"""
    Particle(q, m)

A charged species' identity: signed charge `q` [C] and mass `m` [kg]. Prefer
[`Electron`](@ref), [`Proton`](@ref), [`Ion`](@ref) over raw SI numbers.
"""
struct Particle{Q,M}
    q::Q   # signed charge [C]
    m::M   # mass [kg]
end

Particle(; z=1, A=1, m=nothing) = Particle(z * E_SI, @something(m, A * MP_SI))

Electron() = Particle(-E_SI, ME_SI)
Proton() = Particle(E_SI, MP_SI)


"""
    Species(particle::Particle, vdf; n)

One physical population: identity (`particle`), kinetic model (`vdf`, a normalized shape
with speeds in `v/c`), and number density `n` (SI m⁻³, or Unitful via the extension).
Solver-agnostic: it carries NO normalization (no `B0`, no reference). This is the shared
"species" that every solver agrees on.
"""
struct Species{P<:Particle,N,V}
    particle::P
    n::N
    vdf::V
end
Species(p::Particle, vdf; n) = Species(p, n, vdf)

"""
    Plasma(species::Species...; B0)

The physical system: its `species` and the ambient magnetic field `B0` (SI Tesla, or
Unitful). `B0` is plasma-global — it sets the scale a solver normalizes to, together
with the solver's own choice of reference. Iterable over its species.
"""
struct Plasma{S,B} <: AbstractPlasma
    species::S
    B0::B
end
Plasma(species::Species...; B0) = Plasma(species, B0)

Base.iterate(p::AbstractPlasma, state=1) = iterate(p.species, state)
Base.length(p::AbstractPlasma) = length(p.species)

include("accessor.jl")
end
