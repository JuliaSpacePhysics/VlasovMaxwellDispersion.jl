function residual(prob::DispersionProblem, ω, k = prob.k)
    p, u = _normalized(prob)
    return residual(p.plasma, ω / u.omega, _scale_k(k, u.k, /); closure = prob.closure)
end

include("solver/muller.jl")
include("solver/GRPF.jl")
include("solver/Continuation.jl")
include("solver/AAA.jl")
include("solver/survey.jl")

"""
    solve(prob::DispersionProblem, [alg]; kw...)
    init(prob::DispersionProblem, [alg]; kw...)

Solve `prob` on its normalized form and report in its units.

A `Seed` with a `Wavevector` refines one root (default [`Muller`](@ref)). With any other `k` —
an ordered wavevector list, or a geometry with one swept axis — a `Seed` continues the branch
with [`Continuation`](@ref) (default; another `alg` becomes its `base`), reporting a root at
each `k` given. Both return a [`DispersionSolution`](@ref).

A `Region` finds all root [`DispersionBranch`](@ref)es of the deflated `det(ω̃²𝒟)` as a
[`SurveySolution`](@ref): `alg` ([`AAA`](@ref) default, or [`GRPF`](@ref)) runs at each point of
the geometry's parameter grid, and per-point roots are linked into sheets by [`link`](@ref),
tuned via `linking` options (`gate` defaults to ⅛ of the box diagonal). The ω box is a soft
window tracked `pad` past every edge. Fixed `k` gives single-point branches. `refine` (default
[`Muller`](@ref); `nothing` keeps raw fit/mesh roots) polishes each root and filters out
candidates with no nearby zero of the det.
"""
function CommonSolve.solve(prob::DispersionProblem, alg = _default_alg(prob); kwargs...)
    p, u = _normalized(prob)
    return _restore(_solve(p, alg; kwargs...), prob, u)
end

function CommonSolve.init(prob::DispersionProblem, alg = _default_alg(prob); kwargs...)
    p, u = _normalized(prob)
    return ProblemCache(_init(p, alg; kwargs...), prob, u)
end

struct ProblemCache{C, P, U}
    cache::C
    prob::P
    units::U
end
CommonSolve.solve!(c::ProblemCache) = _restore(solve!(c.cache), c.prob, c.units)
CommonSolve.step!(c::ProblemCache) = (step!(c.cache); c)

_default_alg(::DispersionProblem{<:Seed, <:Wavevector}) = Muller()
_default_alg(::DispersionProblem{<:Seed}) = Continuation()
_default_alg(::DispersionProblem{<:Region}) = AAA()

# Solvers implement `_init` (or `_solve`) on the normalized problem.
_solve(prob, alg; kwargs...) = solve!(_init(prob, alg; kwargs...))
_init(prob::DispersionProblem{<:Seed}, alg; kwargs...) = _init(prob, Continuation(base = alg); kwargs...)

_restore(sol::DispersionSolution, prob, u) =
    DispersionSolution(sol.omega .* u.omega, sol.resid, sol.stats, sol.retcode, prob, sol.alg)
function _restore(sol::SurveySolution, prob, u)
    roots = map(b -> DispersionBranch(b.omega .* u.omega, _scale_k(b.k, u.k), b.resid, b.omega_unit * u.omega), sol.roots)
    return SurveySolution(roots, sol.stats, sol.retcode, prob, sol.alg)
end

function _in_box(region, point = 0)
    ll, ur = region
    return real(ll) ≤ real(point) ≤ real(ur) && imag(ll) ≤ imag(point) ≤ imag(ur)
end
