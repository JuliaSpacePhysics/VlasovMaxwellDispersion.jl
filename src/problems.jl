# SciML-style problem/algorithm/solution layer over the root problem det(𝒟)=0.

abstract type AbstractDispersionProblem end

"""
    DispersionProblem(plasma, omega0, geometry; closure=HarmonicSum())

Refine a seeded mode of `det(𝒟(plasma,ω,k))=0`. A [`Wavenumber`](@ref) `geometry`
selects point refinement. An ordered wavenumber collection, or a geometry selects
[`Continuation`](@ref) along that path.
"""
struct DispersionProblem{P,K,T,C} <: AbstractDispersionProblem
    plasma::P
    omega0::T
    k::K
    closure::C
end
DispersionProblem(plasma, omega0, geometry; closure=HarmonicSum()) =
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
struct GlobalDispersionProblem{P,B,G,C} <: AbstractDispersionProblem
    plasma::P
    region::B
    geometry::G
    closure::C
end
GlobalDispersionProblem(plasma, region, geometry; closure=HarmonicSum()) =
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


module ReturnCode
"""
    ReturnCode.T

Return codes are notes given by the solvers to indicate the state of the solution.

- `Success`   — every requested root converged.
- `Partial`   — a branch was tracked but some `k` are `NaN`.
- `Saturated` — a survey's rational fit stopped short of `tol` (degree cap or
                stagnation), so it may have **missed roots**.
- `MaxIters`  — the polisher hit its iteration cap without reaching tolerance.
- `Failure`   — no root found at all.
"""
@enum T Success Partial Saturated MaxIters Failure
end

"""
    successful_retcode(x)::Bool

True when `x` reports a fully converged solve.
"""
successful_retcode(c::ReturnCode.T) = c === ReturnCode.Success
successful_retcode(sol) = successful_retcode(sol.retcode)

"""
    DispersionSolution

Result of a seeded `solve`. `omega` is a root for point refinement or a vector
of roots for continuation.
`resid` is the scale-invariant [`residual`](@ref) at the root(s),
mirroring the shape of `omega` (`NaN` for non-converged entries).
"""
struct DispersionSolution{T,R,S,Pr,A}
    omega::T
    resid::R
    stats::S
    retcode::ReturnCode.T
    prob::Pr
    alg::A
end


"""
Dispersion branch for parameter sweep. 
For a fixed wavevector, `omega`, `k`, and `resid` are scalars. 
Otherwise they are arrays shaped like the parameter grid.
Missing branch points are `NaN` in `omega` and `resid`.
"""
struct DispersionBranch{W,K,R}
    omega::W
    k::K
    resid::R
end

Base.length(b::DispersionBranch) = length(b.omega)
Base.iterate(b::DispersionBranch, args...) = iterate(b.omega, args...)
Base.getindex(b::DispersionBranch, args...) = getindex(b.omega, args...)

"""
    SurveySolution

A collection of discovered [`DispersionBranch`](@ref)es.
"""
struct SurveySolution{BR,S,Pr,A} <: AbstractVector{BR}
    roots::Vector{BR}
    stats::S
    retcode::ReturnCode.T
    prob::Pr
    alg::A
end

# Prune branches post hoc, e.g. filter(b -> count(isfinite, b.omega) ≥ 5, sol)
Base.filter(pred, s::SurveySolution) = SurveySolution(
    filter(pred, s.roots), s.stats, s.retcode, s.prob, s.alg
)
Base.size(s::SurveySolution) = size(s.roots)
Base.getindex(s::SurveySolution, args...) = s.roots[args...]

_show_stats(io, st::SolveStats) =
    print(io, st.nevals, " evals, ", round(st.time; sigdigits=3), " s")

function Base.show(io::IO, s::SurveySolution)
    print(io, "SurveySolution(retcode=", s.retcode, ", ", length(s.roots), " roots, ")
    _show_stats(io, s.stats)
    return print(io, ")")
end

function Base.show(io::IO, sol::DispersionSolution)
    print(io, "DispersionSolution(retcode=", sol.retcode, ", omega=", sol.omega, ", ")
    _show_stats(io, sol.stats)
    return print(io, ")")
end

Base.show(io::IO, ::MIME"text/plain", sol::SurveySolution) = show(io, sol)
