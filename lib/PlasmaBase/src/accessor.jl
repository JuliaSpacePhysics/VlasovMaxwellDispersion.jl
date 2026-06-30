charge(p::Particle) = p.q
mass(p::Particle) = p.m

frequency(x) = x
magnetic_field(x) = x

particle(s::Species) = s.particle
number_density(s::Species) = s.n
distribution(s::Species) = s.vdf
charge(s::Species) = charge(s.particle)
mass(s::Species) = mass(s.particle)
species(p::Plasma) = p.species
magnetic_field(p::Plasma) = p.B0

"Signed dimensionless gyrofrequency ratio `Ω_s/Ω_ref = (q_s/|q_ref|)(m_ref/m_s)`. B-free."
gyrofrequency_ratio(p::Particle, ref::Particle) = (charge(p) / abs(charge(ref))) * (mass(ref) / mass(p))

"""
    plasma_gyro_ratio(n, m, B0)

Reference-free `ω_ps/Ω_s = c/v_A,s = √(n m/ε₀)/B0` (charge-independent). SI: `n`[m⁻³],
`m`[kg], `B0`[T]; Unitful via the extension.
"""
plasma_gyro_ratio(n, m, B0) = sqrt(n * m / EPS0_SI) / B0
