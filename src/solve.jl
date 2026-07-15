residual(prob::AbstractDispersionProblem, ω, k = prob.k) =
    residual(prob.plasma, ω, k; closure = prob.closure)

include("solver/muller.jl")
include("solver/GRPF.jl")
include("solver/Continuation.jl")
include("solver/AAA.jl")
include("solver/survey.jl")

"""
    solve(prob::DispersionProblem, alg = Muller()) -> DispersionSolution

A [`Wavenumber`](@ref) refines a single seeded root. Any other `k` — an ordered
wavenumber list, or a geometry with one swept axis — continues the branch with
[`Continuation`](@ref), reporting a root at each `k` given.
"""
CommonSolve.solve(prob::DispersionProblem{<:Any, <:Wavenumber}) = CommonSolve.solve(prob, Muller())
CommonSolve.solve(prob::DispersionProblem) = CommonSolve.solve(prob, Continuation())
CommonSolve.init(prob::DispersionProblem, alg; kwargs...) =
    CommonSolve.init(prob, Continuation(base = alg); kwargs...)

"""
    solve(prob::GlobalDispersionProblem, alg=AAA(); refine=Muller(), kw...)::SurveySolution

Find all root [`DispersionBranch`](@ref)es of the deflated `det(ω̃²𝒟)`: `alg`
([`AAA`](@ref) or [`GRPF`](@ref)) runs at each point of the geometry's
parameter grid, and per-point roots are linked into sheets by
[`link`](@ref), tuned via `linking` options (`gate` defaults to ⅛ of the box diagonal).
The ω box is a soft window tracked `pad` past every edge. 
Fixed `k` gives single-point branches.

`refine` (default [`Muller`](@ref); `nothing` keeps raw fit/mesh roots)
polishes each root and filters out candidates with no nearby zero of the det.
"""
CommonSolve.solve(prob::GlobalDispersionProblem; kwargs...) =
    CommonSolve.solve(prob, AAA(); kwargs...)

function _in_box(region, point = 0)
    ll, ur = region
    return real(ll) ≤ real(point) ≤ real(ur) && imag(ll) ≤ imag(point) ≤ imag(ur)
end
