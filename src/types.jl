# Orbit-integral evaluator (derivation.md §3): the `closure` of `AdaptiveEval` and of the direct `contribution`.
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

@inline vec3(k::Wavevector) = SVector(k.k_perp, zero(k.k_perp), k.k_para)

_complex_nan(x) = oftype(complex(x), complex(NaN, NaN))

regime(d::NormalizedPopulation) = regime(d.vdf)

# This solver's form of a population: `vdf` converted by `solver_vdf`, speeds in `c`. Idempotent,
# so every population-level entry (`contribution`, `plan_contribution`) converts on the way in.
_solver_population(s) = NormalizedPopulation(s.Omega, s.Pi2, solver_vdf(s.vdf, s.c), 1)

prepare(s::NormalizedPopulation; kw...) =
    (s = _solver_population(s); NormalizedPopulation(s.Omega, s.Pi2, prepare(s.vdf; kw...), 1))
prepare(p::NormalizedPlasma; kw...) = NormalizedPlasma(map(s -> prepare(s; kw...), p.populations), p.scales)
prepare(p::Plasma; kw...) = prepare(NormalizedPlasma(p); kw...)

"""VDF spec plus its `precompute`d (ω,k)-independent constants."""
struct PreparedVDF{V,C} <: SolverVDF
    vdf::V
    cache::C
end

regime(c::PreparedVDF) = regime(c.vdf)
trusted(c::PreparedVDF, s, ω, k) = trusted(c.vdf, s, ω, k)
parallel_even(c::PreparedVDF) = parallel_even(c.vdf)
