# (B0, n, q, m, T) → the dimensionless scales a solver normalizes against.
# `B0` is SI here: the `Plasma` constructor strips units once, at the boundary.

"""
    Scales

A plasma's dimensionless scales, as built by [`scales`](@ref). Property access forwards to
the underlying `NamedTuple`; the type exists so a VDF can complete itself by dispatch,
`MyVDF(sc::Scales)`.
"""
struct Scales{NT<:NamedTuple}
    nt::NT
end

Base.getproperty(sc::Scales, k::Symbol) = getproperty(getfield(sc, :nt), k)
Base.propertynames(sc::Scales) = propertynames(getfield(sc, :nt))
Base.NamedTuple(sc::Scales) = getfield(sc, :nt)
Base.show(io::IO, sc::Scales) = print(io, "Scales", NamedTuple(sc))

"""
    scales(plasma)    -> Scales
    scales(plasma, i) -> Scales

The plasma's characteristic scales, dimensionless — speeds in `c`, lengths in `c/Ω_ref`.

- `vA` — Alfvén speed, from the total mass density
- `Omega_ref`, `B0` — reference frequency [rad/s] and field [T] fixing the normalization
- `species` — per species `(; wc, wp, vth_perp, vth_para, n, d, rho)`, with normalized
  gyro- and plasma frequency `wc = Ω_s/Ω_ref` (signed) and `wp = ω_ps/Ω_ref`, inertial
  length `d = c/ω_ps = 1/wp` and gyroradius `rho = vth_perp/|Ω_s|`. Both axes are named
  explicitly; anisotropy pairs are ordered `(⊥, ∥)` — see [`perp_para`](@ref)

`scales(plasma, i)` merges the `i`-th species' entry with the plasma-wide ones.
"""
scales(p::Plasma) = _scales(species(p), magnetic_field(p), reference(p))
scales(sc::Scales) = sc
scales(p, i::Integer) = _species_view(scales(p), i)

_species_view(sc::Scales, i) = Scales(_species_view(NamedTuple(sc), i))
_species_view(sc::NamedTuple, i) = merge(Base.structdiff(sc, (; species=nothing)), sc.species[i])

function _scales(sps, B0, ref)
    sc = map(s -> _species_scales(s, B0, ref), Tuple(sps))
    Omega_ref = ref isa Particle ? abs(charge(ref)) * B0 / mass(ref) : frequency(ref)
    vA = 1 / sqrt(sum(s -> (s.wp / s.wc)^2, sc))
    return Scales((; vA, B0, Omega_ref, species=sc))
end

function _species_scales(s::Species, B0, ref)
    p = particle(s)
    wc = ref isa Particle ? gyrofrequency_ratio(p, ref) : (charge(p) * B0 / mass(p)) / frequency(ref)
    wp = plasma_gyro_ratio(number_density(s), mass(p), B0) * abs(wc)
    vth_perp, vth_para = thermal_speeds(s, B0)
    return (;
        wc, wp, vth_perp, vth_para, n=number_density(s), d=1 / wp,
        rho=isnothing(vth_perp) ? nothing : vth_perp / abs(wc),
    )
end

"""
    thermal_speeds(s::Species, B0) -> (vth_perp, vth_para)

The species' thermal speeds in `v/c`, from whichever of `T`/`beta`/`vth` it carries;
`(nothing, nothing)` when it declares no thermal scale.
"""
function thermal_speeds(s::Species, B0)
    th = s.thermal
    haskey(th, :vth) && return map(velocity, perp_para(th.vth))
    haskey(th, :T) && return map(T -> sqrt(2 * energy(T) / mass(s)) / C_SI, perp_para(th.T))
    haskey(th, :beta) || return (nothing, nothing)
    # β = 2μ₀nT/B² ⇒ vth = √(2T/m) = √β·(B/√(μ₀ n m)), the species' own Alfvén speed
    vA_s = magnetic_field(B0) / sqrt(MU0_SI * mass_density(s)) / C_SI
    return map(β -> sqrt(β) * vA_s, perp_para(th.beta))
end

"""
    perp_para(x) -> (perp, para)

Widen an anisotropy parameter to its `(⊥, ∥)` pair: a scalar is isotropic, a
`Tuple`/`AbstractVector` is positional, a `NamedTuple` is order-free.
"""
perp_para(x) = (x, x)
perp_para(x::Union{Tuple,AbstractVector}) = (x[1], x[2])
perp_para(x::NamedTuple) = (x.perp, x.para)
