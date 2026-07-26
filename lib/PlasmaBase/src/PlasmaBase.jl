"""
    PlasmaBase

Thin, solver-agnostic vocabulary for describing a plasma *physically*.
It holds the PROBLEM description (identity, populations, the system and its field) plus an accessor interface.
"""
module PlasmaBase

export AbstractVDF, Particle, Electron, Proton, Species, Plasma
export charge, mass, particle, number_density, distribution, species, magnetic_field
export gyrofrequency_ratio, plasma_gyro_ratio
export Scales, scales

"""Supertype of every velocity distribution function."""
abstract type AbstractVDF end
abstract type AbstractPlasma end

# --- SI constants (plain Float64) ---
const C_SI = 2.99792458e8       # m/s
const E_SI = 1.602176634e-19    # C
const ME_SI = 9.1093837015e-31   # kg
const MP_SI = 1.67262192369e-27  # kg
const EPS0_SI = 8.8541878128e-12   # F/m
const MU0_SI = 1 / (EPS0_SI * C_SI^2)  # H/m
const KB_SI = 1.380649e-23       # J/K

"""
    Particle(q, m)

A charged species' identity: signed charge `q` [C] and mass `m` [kg]. 
"""
struct Particle{Q,M}
    q::Q   # signed charge [C]
    m::M   # mass [kg]
end

"""
    Particle(; z=1, A=1, m=nothing)

Particle with charge number `z` and mass number `A`.
"""
Particle(; z=1, A=1, m=nothing) = Particle(z * E_SI, @something(m, A * MP_SI))

Electron() = Particle(-E_SI, ME_SI)
Proton() = Particle(E_SI, MP_SI)


"""
    Species(particle, vdf=nothing; n, T=nothing, beta=nothing, vth=nothing)

The thermal input is at most one of `T` (temperature), `beta` (plasma beta), or `vth` (thermal speed).
Scalar for isotropic, tuple for anisotropic (parallel, perpendicular).
"""
struct Species{P,N,V,K}
    particle::P
    n::N
    vdf::V
    thermal::K   # whichever of (T, beta, vth) the user supplied
end

const THERMAL_KEYS = (:T, :beta, :vth)

function Species(p, vdf=nothing; n, kw...)
    th = NamedTuple(kw)
    bad = filter(∉(THERMAL_KEYS), keys(th))
    isempty(bad) && length(th) ≤ 1 ||
        throw(ArgumentError("Species accepts `n` and at most one of $THERMAL_KEYS, got $(keys(th))"))
    return Species(p, number_density(n), vdf, th)
end

"""
    Plasma(species::Species...; B0, ref=first species' particle)

The physical system: its `species`, the ambient magnetic field `B0` (SI Tesla, or
Unitful), and the reference `ref` that fixes the normalization — a `Particle`, whose
gyrofrequency becomes `Ω_ref`, or any frequency in rad/s. `ref` need not be one of the
species. Iterable over its species.
"""
struct Plasma{S,B,R} <: AbstractPlasma
    species::S
    B0::B
    ref::R
end
Plasma(species::Species...; B0, ref=nothing) =
    Plasma(species, magnetic_field(B0), @something(ref, particle(first(species))))

Base.iterate(p::AbstractPlasma, state=1) = iterate(p.species, state)
Base.length(p::AbstractPlasma) = length(p.species)

include("accessor.jl")
include("scales.jl")
end
