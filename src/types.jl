abstract type AbstractDispersionProblem end

# Orbit-integral evaluator (derivation.md §3), passed as the `closure` keyword of contribution/solve.
abstract type IntegralClosure end
struct HarmonicSum <: IntegralClosure end
struct Newberger <: IntegralClosure end

abstract type Regime end
struct NonRelativistic <: Regime end
struct Relativistic <: Regime end

# ── Swept manifold in k-space
abstract type ParameterGeometry end

"""
    Wavenumber(kperp, kz)
    Wavenumber(; kz, kperp=zero(kz))

Dimensionless wavevector `k c / Omega_ref`, `(kperp, kz)`.
"""
struct Wavenumber{T}
    kperp::T
    kz::T
end
Wavenumber(kperp, kz) = Wavenumber(promote(kperp, kz)...)
Wavenumber(; kz, kperp = zero(kz)) = Wavenumber(kperp, kz)

Base.eltype(::Type{Wavenumber{T}}) where {T} = T
Base.convert(::Type{Wavenumber{T}}, k::Wavenumber) where {T} = Wavenumber{T}(k.kperp, k.kz)
Base.broadcastable(k::Wavenumber) = Ref(k)

@inline para(k::Wavenumber) = k.kz
@inline perp(k::Wavenumber) = k.kperp
@inline Base.abs2(k::Wavenumber) = k.kz^2 + k.kperp^2
@inline Base.angle(k::Wavenumber) = atan(k.kperp, k.kz)  # propagation angle to B0
@inline vec3(k::Wavenumber) = SVector(k.kperp, zero(k.kperp), k.kz)

"""
    AngleSweep(k, theta)   /   AngleSweep(; k, theta)

Polar k-space geometry: `|k|` and the angle `θ` to `B₀` so that `𝐤 = (k sinθ, k cosθ)`.
Each axis is a fixed scalar or swept — a vector of sample values, or a `(lo, hi)` tuple.
"""
struct AngleSweep{K, T} <: ParameterGeometry
    k::K
    theta::T
end
AngleSweep(; k, theta) = AngleSweep(k, theta)

"""CartesianSweep(kx, kz)   /   CartesianSweep(; kx=false, kz)"""
struct CartesianSweep{P, Z} <: ParameterGeometry
    kx::P
    kz::Z
end
# kx=false: Bool is the weakest promoting Number, so kz's type wins.
CartesianSweep(; kx = false, kz) = CartesianSweep(kx, kz)

_as_tuple(x::T) where {T} = Tuple{fieldtypes(T)...}(getfield(x, i) for i in 1:fieldcount(T))

_isswept(x) = x isa Union{Tuple, AbstractVector}
_wavenumber(::AngleSweep, k, θ) = Wavenumber(k .* sincos(θ)...)
_wavenumber(::CartesianSweep, kx, kz) = Wavenumber(kx, kz)
_wavenumber(::Wavenumber, kperp, kz) = Wavenumber(kperp, kz)

# fixed axis emits its scalar; swept axis takes the next param, returning the rest.
_pick(axis::Number, ps) = (axis, ps)
_pick(axis, ps) = (first(ps), Base.tail(ps))

function wavefun(g)
    a, b = _as_tuple(g)
    return function (ps...)
        va, rest = _pick(a, ps)
        vb, _ = _pick(b, rest)
        return _wavenumber(g, va, vb)
    end
end

_paramgrid(axis::AbstractVector) = axis
_paramgrid(axis::Tuple) = range(axis...; length = 61)
paramgrids(g) = map(_paramgrid, filter(_isswept, _as_tuple(g)))

# Real float type carried by a value/type — the element-type anchor so user
# number types (Float32, BigFloat, …) flow into solution containers.
_realtype(x) = float(real(eltype(x)))
_realtype(x::Tuple) = promote_type(map(_realtype, x)...)
_realtype(g::ParameterGeometry) = _realtype(_as_tuple(g))

_complex_nan(x) = oftype(complex(x), complex(NaN, NaN))

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


prepare(s::NormalizedSpecies, args...; kw...) =
    NormalizedSpecies(s.Omega, s.Pi2, prepare(s.vdf, args...; kw...))
prepare(p::NormalizedPlasma, args...; kw...) =
    NormalizedPlasma(map(s -> prepare(s, args...; kw...), p.species))

"""VDF spec plus its `precompute`d (ω,k)-independent constants."""
struct PreparedVDF{V, C} <: AbstractVDF
    vdf::V
    cache::C
end

regime(c::PreparedVDF) = regime(c.vdf)
trusted(c::PreparedVDF, s, ω, k) = trusted(c.vdf, s, ω, k)
parallel_even(c::PreparedVDF) = parallel_even(c.vdf)
