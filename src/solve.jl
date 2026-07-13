dispersion_function(prob::AbstractDispersionProblem) =
    ω -> det(dispersion_tensor(prob.plasma, ω, prob.k; closure = prob.closure))

# det(ω̃²𝒟)=ω̃⁶·det𝒟 is pole-free at ω=0. The genuine light-term pole there
# would otherwise cancel nearby roots' winding
deflated_dispersion_function(prob::AbstractDispersionProblem) =
    ω -> det(wave_dispersion_tensor(prob.plasma, ω, prob.k; closure = prob.closure))

@inline Base.getproperty(prob::DispersionProblem, s::Symbol) =
    s === :f ? _scaled_dispersion_function(prob) : getfield(prob, s)
@inline Base.getproperty(prob::GlobalDispersionProblem, s::Symbol) =
    s === :f ? deflated_dispersion_function(prob) : getfield(prob, s)
Base.propertynames(prob::AbstractDispersionProblem) =
    (fieldnames(typeof(prob))..., :f)

_hadamard(D) = prod(norm(D[i, :]) for i in 1:3)

residual(prob::AbstractDispersionProblem, ω, k = prob.k) =
    residual(prob.plasma, ω, k; closure = prob.closure)

function _scaled_dispersion_function(prob::DispersionProblem)
    f = dispersion_function(prob)
    s = _hadamard(dispersion_tensor(prob.plasma, prob.omega0, prob.k; closure = prob.closure))
    return isfinite(s) && s > 0 ? (ω -> f(ω) / s) : f
end

include("solver/muller.jl")
include("solver/GRPF.jl")
include("solver/ArcLength.jl")
include("solver/AAA.jl")
include("solver/survey.jl")

"""
    solve(prob::DispersionProblem, alg = Muller()) -> DispersionSolution
"""
CommonSolve.solve(prob::DispersionProblem{<:Any, <:Wavenumber}) = CommonSolve.solve(prob, Muller())
CommonSolve.solve(prob::DispersionProblem) = CommonSolve.solve(prob, ArcLength())
CommonSolve.init(prob::DispersionProblem, alg; kwargs...) =
    CommonSolve.init(prob, ArcLength(base = alg); kwargs...)

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
