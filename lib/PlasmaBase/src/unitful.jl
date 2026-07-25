# Unitful input adapters for the physical vocabulary
#
# Plain numbers stay valid (SI base units, eV for temperature, v/c for speed); Unitful
# lets you write 5u"cm^-3", 5u"nT", 10u"eV", 400u"km/s". Everything converts to SI Float64.

using Unitful
using Unitful: AbstractQuantity, ustrip, uconvert, dimension, @u_str

_kg(m) = m isa AbstractQuantity ? ustrip(u"kg", m) : m

# temperature units carry no energy dimension until multiplied by k_B
_joule(x::AbstractQuantity) =
    dimension(x) == dimension(u"K") ? ustrip(u"J", uconvert(u"J", x * Unitful.k)) : ustrip(u"J", uconvert(u"J", x))

PB.magnetic_field(B::AbstractQuantity) = ustrip(u"T", B)
PB.frequency(ω::AbstractQuantity) = ustrip(u"s^-1", uconvert(u"s^-1", ω))
PB.number_density(n::AbstractQuantity) = ustrip(u"m^-3", n)
PB.energy(T::AbstractQuantity) = _joule(T)
PB.velocity(v::AbstractQuantity) = ustrip(u"m/s", uconvert(u"m/s", v)) / PB.C_SI

PB.plasma_gyro_ratio(n::AbstractQuantity, m, B::AbstractQuantity) =
    PB.plasma_gyro_ratio(ustrip(u"m^-3", n), _kg(m), ustrip(u"T", B))

PB.Particle(q::AbstractQuantity, m::AbstractQuantity) =
    PB.Particle(ustrip(u"C", q), ustrip(u"kg", m))
