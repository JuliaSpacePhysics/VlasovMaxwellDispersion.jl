# Bridge: normalize a PHYSICAL `PlasmaBase.Plasma` into this solver's dimensionless rep.
# map (B0, n, q, m) → (Omega, Pi2)

# `ref` a frequency in rad/s ⇒ Omega = q B0/m / Ω_ref (carries B0).
# Ω_ref is taken POSITIVE (via `abs(q_ref)` in gyrofrequency_ratio) so sign(Omega)=sign(q_s).
_normalized_omega(p::Particle, ref::Particle, _B0) = gyrofrequency_ratio(p, ref)
_normalized_omega(p::Particle, Ωref, B0) = (charge(p) * magnetic_field(B0) / mass(p)) / frequency(Ωref)

"""
    NormalizedSpecies(s::Species, B0, ref=particle(s))

Normalize one physical `Species` against the reference `ref` and field `B0`:
`Omega = Ω_s/Ω_ref` and `Pi2 = (ω_ps/Ω_ref)² = (wpwc·Omega)²`, `wpwc = √(n m/ε₀)/B0`.
`ref` is a [`Particle`](@ref) (Ω_ref = its gyrofrequency, Omega is B-free) or any
reference frequency.
"""
function NormalizedSpecies(s::Species, B0::Number, ref = particle(s))
    p = particle(s)
    Ω = _normalized_omega(p, ref, B0)
    wpwc = plasma_gyro_ratio(number_density(s), mass(p), B0)
    return NormalizedSpecies(Ω, (wpwc * abs(Ω))^2, distribution(s))
end


"""
    NormalizedPlasma(plasma::Plasma; ref=first particle)

Normalize a physical plasma into this solver's dimensionless [`NormalizedPlasma`](@ref).
`B0` comes from the plasma; the reference `ref` (a [`Particle`](@ref) or a reference
frequency) sets `Ω_ref` and defaults to the first species' particle.
"""
function NormalizedPlasma(plasma::Plasma; ref = particle(first(plasma)))
    B0 = magnetic_field(plasma)
    return NormalizedPlasma(map(s -> NormalizedSpecies(s, B0, ref), plasma))
end

prepare(p::Plasma, args...; kw...) = prepare(NormalizedPlasma(p), args...; kw...)
