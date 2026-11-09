# Orbit-integral evaluator (derivation.md §3), passed as the `closure` keyword of contribution/solve.
abstract type IntegralClosure end
struct HarmonicSum <: IntegralClosure end
struct Newberger <: IntegralClosure end

abstract type Regime end
struct NonRelativistic <: Regime end
struct Relativistic <: Regime end

"""
    SolverVDF

Supertype of VDFs in the solver's own representation, momenta in units of `mc`.
"""
abstract type SolverVDF end

"""
    Seed(omega0, at=nothing)

A single root near `omega0` where optional `at` anchors it to a wavevector.
"""
struct Seed{W,K}
    omega0::W
    at::K
end
Seed(omega0) = Seed(omega0, nothing)
Base.getindex(s::Seed, args...) = Base.getindex(s.omega0, args...)

struct Region{B}
    box::B
end

Base.iterate(r::Region, state=1) = iterate(r.box, state)
Base.getindex(r::Region, i) = r.box[i]

@inline vec3(k::Wavevector) = SVector(k.k_perp, zero(k.k_perp), k.k_para)

_as_tuple(x::T) where {T} = Tuple{fieldtypes(T)...}(getfield(x, i) for i in 1:fieldcount(T))

_isswept(::AbstractVector) = true
_isswept(::Number) = false
_isswept(x) = throw(ArgumentError("sweep axis must be a number or a vector of values, got $(typeof(x))"))
_wavevector(::AngleSweep, k, θ) = Wavevector(k .* sincos(θ)...)
_wavevector(::CartesianSweep, kperp, kpara) = Wavevector(kperp, kpara)
_wavevector(::Wavevector, kperp, kpara) = Wavevector(kperp, kpara)

# fixed axis emits its scalar; swept axis takes the next param, returning the rest.
_pick(axis::Number, ps) = (axis, ps)
_pick(::AbstractVector, ps) = (first(ps), Base.tail(ps))

function wavefun(g)
    a, b = _as_tuple(g)
    return function (ps...)
        va, rest = _pick(a, ps)
        vb, _ = _pick(b, rest)
        return _wavevector(g, va, vb)
    end
end

paramgrids(g) = filter(_isswept, _as_tuple(g))

# Real float type carried by a value/type — the element-type anchor so user
# number types (Float32, BigFloat, …) flow into solution containers.
_realtype(x) = float(real(eltype(x)))
_realtype(x::Tuple) = promote_type(map(_realtype, x)...)
_realtype(g::ParameterGeometry) = _realtype(_as_tuple(g))

_complex_nan(x) = oftype(complex(x), complex(NaN, NaN))

"""
    NormalizedPopulation(Omega, Pi2, vdf)

Solver's dimensionless per-population representation. `Omega = Ω_s/Ω_ref`; `Pi2 = (ω_ps/Ω_ref)^2`, negative for a
negative-density population.
A PlasmaPopulations `vdf` takes speeds in units of `c` here and is converted by [`solver_vdf`](@ref).
"""
struct NormalizedPopulation{T,V}
    Omega::T
    Pi2::T
    vdf::V
    global _normalized_population(Omega::T, Pi2::T, vdf::V) where {T,V} = new{T,V}(Omega, Pi2, vdf)
end

NormalizedPopulation(Omega, Pi2, vdf) = _normalized_population(Omega, Pi2, vdf, 1)
_normalized_population(Omega, Pi2, vdf, vc) =
    _normalized_population(promote(Omega, Pi2)..., solver_vdf(vdf, vc))

regime(d::NormalizedPopulation) = regime(d.vdf)

"""
    NormalizedPlasma(populations...)
    NormalizedPlasma(plasma::Plasma; ref)

Solver's dimensionless container: a bag of [`NormalizedPopulation`](@ref), buildable from a
physical `Plasma` with frequencies normalized to `Ω_ref`, see [`scales`](@ref).
"""
struct NormalizedPlasma{S,C}
    populations::S
    scales::C
end
NormalizedPlasma(populations) = NormalizedPlasma(populations, nothing)
NormalizedPlasma(populations::NormalizedPopulation...) = NormalizedPlasma(Tuple(populations), nothing)
NormalizedPlasma(p::NormalizedPlasma) = p

prepare(s::NormalizedPopulation, args...; kw...) =
    _normalized_population(s.Omega, s.Pi2, prepare(s.vdf, args...; kw...))
prepare(p::NormalizedPlasma, args...; kw...) =
    NormalizedPlasma(map(s -> prepare(s, args...; kw...), p.populations), p.scales)

"""VDF spec plus its `precompute`d (ω,k)-independent constants."""
struct PreparedVDF{V,C} <: SolverVDF
    vdf::V
    cache::C
end

regime(c::PreparedVDF) = regime(c.vdf)
trusted(c::PreparedVDF, s, ω, k) = trusted(c.vdf, s, ω, k)
parallel_even(c::PreparedVDF) = parallel_even(c.vdf)
