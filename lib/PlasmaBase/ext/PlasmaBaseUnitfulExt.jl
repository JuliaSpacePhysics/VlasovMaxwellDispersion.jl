module PlasmaBaseUnitfulExt

# Unitful input adapter for the physical helpers: dimensional safety where it pays (the
# n×B group and the thermal speed). Plain numbers remain valid as SI base units; Unitful
# lets you write 5u"cm^-3", 5u"nT", 10u"eV". Both convert to SI Float64 and call core.
import PlasmaBase as PB
using Unitful
using Unitful: AbstractQuantity

_kg(m) = m isa AbstractQuantity ? ustrip(u"kg", m) : m

PB.magnetic_field(B::AbstractQuantity) = ustrip(u"T", B)
PB.frequency(ω::AbstractQuantity) = ustrip(u"s^-1", uconvert(u"s^-1", ω))

PB.plasma_gyro_ratio(n::AbstractQuantity, m, B::AbstractQuantity) =
    PB.plasma_gyro_ratio(ustrip(u"m^-3", n), _kg(m), ustrip(u"T", B))

PB.Particle(q::AbstractQuantity, m::AbstractQuantity) =
    PB.Particle(ustrip(u"C", q), ustrip(u"kg", m))

end
