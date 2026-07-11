# SciML-style problem/algorithm/solution layer over the root problem det(𝒟)=0.

abstract type AbstractDispersionProblem end

"""
    DispersionProblem(plasma, omega0, geometry; closure=HarmonicSum())

Refine a seeded mode of `det(𝒟(plasma,ω,k))=0`. A [`Wavenumber`](@ref) `geometry`
selects point refinement; an ordered wavenumber collection selects continuation.
"""
struct DispersionProblem{P, K, T, C} <: AbstractDispersionProblem
    plasma::P
    omega0::T
    k::K
    closure::C
end
DispersionProblem(plasma, omega0, geometry; closure = HarmonicSum()) =
    DispersionProblem(plasma, omega0, geometry, closure)


"""
    GlobalDispersionProblem(plasma, region, geometry; closure=HarmonicSum())

Find *all* root branches of `det(𝒟)` in the complex ω box
`region=(lowerleft, upperright)` as the k-space `geometry` sweeps its parameters. 
A branch is an `m`-manifold `ω*(p)` where `m` is the geometry's swept dimension: 
a point (`m=0`), curve (`m=1`), surface (`m=2`), ….

`geometry` can be:
- a [`Wavenumber`](@ref) — fixed k, `m=0`;
- an [`AngleSweep`](@ref) `(|k|, θ)` or [`CartesianSweep`](@ref) `(k⊥, k∥)`
"""
struct GlobalDispersionProblem{P, B, G, C} <: AbstractDispersionProblem
    plasma::P
    region::B
    geometry::G
    closure::C
end
GlobalDispersionProblem(plasma, region, geometry; closure = HarmonicSum()) =
    GlobalDispersionProblem(plasma, region, geometry, closure)

prepare(prob::DispersionProblem; kw...) =
    DispersionProblem(prepare(prob.plasma, prob.closure; kw...), prob.omega0, prob.k, prob.closure)
prepare(prob::GlobalDispersionProblem; kw...) =
    GlobalDispersionProblem(prepare(prob.plasma, prob.closure; kw...), prob.region, prob.geometry, prob.closure)

_realtype(p::GlobalDispersionProblem) =
    promote_type(_realtype(p.region), _realtype(p.geometry))


struct SolveStats{T}
    nevals::Int # number of evaluations
    time::T # solving time in seconds
end

"""
    DispersionSolution

Result of a seeded `solve`. `omega` is a root for point refinement or a vector
of roots for continuation. `retcode` is
`:Success`, `:Failure`, or `:Partial` (branch with some non-converged `k`).
`resid` is the scale-invariant [`residual`](@ref) at the root(s),
mirroring the shape of `omega` (`NaN` for non-converged entries).
"""
struct DispersionSolution{T, R, S, Pr, A}
    omega::T
    resid::R
    stats::S
    retcode::Symbol
    prob::Pr
    alg::A
end


"""
    DispersionBranch

One surveyed branch: mode `omega` at wavevector `k` with residual `resid`.
"""
struct DispersionBranch{W, K, R}
    omega::W
    k::K
    resid::R
end

Base.getindex(b::DispersionBranch, args...) = getindex(b.omega, args...)

"""
    SurveySolution

A collection of discovered [`DispersionBranch`](@ref)es.
`retcode` is `:Success`, `:Failure` (no root branch found),
or `:Partial` (a fit saturated — structure may exceed one rational fit).
"""
struct SurveySolution{BR, Pr, A}
    roots::Vector{BR}
    stats::SolveStats
    retcode::Symbol
    prob::Pr
    alg::A
end

# Prune branches post hoc, e.g. filter(b -> length(b.omega) ≥ 5, sol)
Base.filter(pred, s::SurveySolution) = SurveySolution(
    filter(pred, s.roots), s.stats, s.retcode, s.prob, s.alg
)

_show_stats(io, st::SolveStats) =
    print(io, st.nevals, " evals, ", round(st.time; sigdigits = 3), " s")

function Base.show(io::IO, s::SurveySolution)
    print(io, "SurveySolution(retcode=:", s.retcode, ", ", length(s.roots), " roots, ")
    _show_stats(io, s.stats)
    return print(io, ")")
end

function Base.show(io::IO, sol::DispersionSolution)
    print(io, "DispersionSolution(retcode=:", sol.retcode, ", omega=", sol.omega, ", ")
    _show_stats(io, sol.stats)
    return print(io, ")")
end
