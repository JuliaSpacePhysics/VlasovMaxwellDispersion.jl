# SciML-style problem/algorithm/solution layer over the root problem det(𝒟)=0.

abstract type AbstractDispersionProblem end

"""
    DispersionProblem(plasma, target, k; closure=HarmonicSum(), mode=:det)

Chase zeros of `det(𝒟(plasma,ω,k))=0`. `target` selects *which* roots:

- [`Seed`](@ref)`(ω, at=nothing)` / a number `ω` — refine/track one branch
- `Region((lowerleft, upperright))` / a tuple — survey every branch in the ω box.

`k` sets the wavenumber domain: a [`Wavenumber`](@ref) (point), an ordered wavenumber
collection (path), or a parameter sweep ([`AngleSweep`](@ref) `(|k|, θ)` 
or [`CartesianSweep`](@ref) `(k⊥, k∥)`).

`mode` selects the [`TensorReduction`](@ref) whose zeros are chased (`:det` default).

`solve` returns [`DispersionSolution`](@ref) for `Seed`, or [`SurveySolution`](@ref)
of all [`DispersionBranch`](@ref)es for `Region` — each an `m`-manifold `ω*(p)` over
the `m` swept axes of `k` (point/curve/surface for `m` = 0/1/2).
"""
struct DispersionProblem{T,K,P,C,M} <: AbstractDispersionProblem
    plasma::P
    target::T
    k::K
    closure::C
    mode::M
end
DispersionProblem(plasma, target, k; closure=HarmonicSum(), mode=FullDet()) =
    DispersionProblem(plasma, _target(target), k, closure, TensorReduction(mode))

function Base.getproperty(prob::DispersionProblem, sym::Symbol)
    sym === :omega0 && return prob.target.omega0
    sym === :f && return DispersionFunction(prob)
    return getfield(prob, sym)
end

_target(x) = x
_target(omega::Number) = Seed(omega)
_target(region::Tuple) = Region(region)

const GlobalDispersionProblem = DispersionProblem

prepare(prob::DispersionProblem; kw...) = DispersionProblem(
    prepare(prob.plasma, prob.closure; kw...), prob.target, prob.k, prob.closure, prob.mode)

_realtype(p::DispersionProblem{<:Region}) =
    promote_type(_realtype(p.target.box), _realtype(p.k))


struct SolveStats{T}
    nevals::Int # number of evaluations
    time::T # solving time in seconds
end

Base.:+(a::SolveStats, b::SolveStats) = SolveStats(a.nevals + b.nevals, a.time + b.time)


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
    SurveySolution

A vector of discovered [`DispersionBranch`](@ref)es: `sol[i]` is the `i`-th
branch, `filter(isgrowing, sol)` keeps the unstable ones. Each branch spans the full
`k` grid (missing points are `NaN`).
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
