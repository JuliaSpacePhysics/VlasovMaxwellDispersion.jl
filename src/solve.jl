dispersion_function(prob::AbstractDispersionProblem) =
    ω -> det(dispersion_tensor(prob.plasma, ω, prob.k; closure = prob.closure))

# det(ω̃²𝒟)=ω̃⁶·det𝒟 is pole-free at ω=0. The genuine light-term pole there
# would otherwise cancel nearby roots' winding
deflated_dispersion_function(prob::AbstractDispersionProblem) =
    ω -> det(wave_dispersion_tensor(prob.plasma, ω, prob.k; closure = prob.closure))

@inline Base.getproperty(prob::LocalDispersionProblem, s::Symbol) =
    s === :f ? _scaled_dispersion_function(prob) : getfield(prob, s)
@inline Base.getproperty(prob::GlobalDispersionProblem, s::Symbol) =
    s === :f ? deflated_dispersion_function(prob) : getfield(prob, s)
Base.propertynames(prob::AbstractDispersionProblem) =
    (fieldnames(typeof(prob))..., :f)

_hadamard(D) = prod(norm(D[i, :]) for i in 1:3)

"""
    residual(plasma, ω, k; closure=HarmonicSum())
    residual(prob, ω)

Scale-invariant residual `|det 𝒟(ω)| / ∏ᵢ‖𝒟ᵢ,:‖ ∈ [0,1]`; ~machine epsilon at a
genuine root regardless of the tensor's magnitude. `NaN` for non-finite `ω`.

Raw |det 𝒟| floors at ~‖𝒟‖³ε from cancellation.
"""
function residual(plasma, ω, k; closure = HarmonicSum())
    isfinite(ω) || return NaN
    D = dispersion_tensor(plasma, ω, k; closure)
    return abs(det(D)) / _hadamard(D)
end
residual(prob::AbstractDispersionProblem, ω, k = prob.k) =
    residual(prob.plasma, ω, k; closure = prob.closure)

function _scaled_dispersion_function(prob::LocalDispersionProblem)
    f = dispersion_function(prob)
    s = _hadamard(dispersion_tensor(prob.plasma, prob.omega0, prob.k; closure = prob.closure))
    return isfinite(s) && s > 0 ? (ω -> f(ω) / s) : f
end

include("solver/muller.jl")
include("solver/GRPF.jl")
include("solver/ArcLength.jl")

"""
    solve(prob::LocalDispersionProblem, alg = Muller()) -> DispersionSolution
"""
CommonSolve.solve(prob::LocalDispersionProblem) = CommonSolve.solve(prob, Muller())


"""
    wave_dispersion_tensor(plasma, ω, k::Wavenumber; closure=HarmonicSum())

Deflated form `ω̃²·𝒟 = ω̃²ε + (k̃k̃ᵀ − k̃²I)`, built as `ω̃²I + ω̃²χ + curlcurl` so 
the light-term `curlcurl/ω̃²` pole for original `det𝒟` and any `χ` pole at `ω=0` 
(cold `ε`'s `1/ω²`, `1/ω` terms) cancel analytically.

Otherwise its winding partially cancels nearby roots, causing GRPF to miss them
and report a spurious net pole.
"""
function wave_dispersion_tensor(plasma, ω, k; kwargs...)
    ω2χ = _guarded_sum(s -> scaled_contribution(s, ω, k; kwargs...), plasma)
    return ω^2 * I + ω2χ + _curlcurl(k)
end

"""
    solve(prob::GlobalDispersionProblem; refine=Muller()) -> SurveySolution

[`SurveySolution`](@ref) contains a list of [`DispersionBranch`](@ref)es (root and pole).

At fixed `k` (`m=0`) each root/pole is a single-point branch.

Global survey (e.g. GRPF) locates roots only to mesh accuracy. 
`refine` method (default: `Muller()`) polishes each root to convergence.
Pass `refine=nothing` to keep the raw mesh roots.
"""
CommonSolve.solve(prob::GlobalDispersionProblem; kwargs...) =
    CommonSolve.solve(prob, GRPF(); kwargs...)
# Fixed-k roots/poles as SurveySolution
function _fixedk_survey(prob, alg, k, roots, poles, nevals, retcode; refine = Muller())
    isnothing(refine) || (roots = _refine_roots(prob, k, roots, refine))
    roots = [DispersionBranch(ω, k, residual(prob, ω, k)) for ω in roots]
    poles = [DispersionBranch(ω, k, nothing) for ω in poles]
    return SurveySolution(roots, poles, nevals, retcode, prob, alg)
end

# Keep mesh value on divergence and drop polished duplicates
function _refine_roots(prob, k, roots, alg)
    polished = eltype(roots)[]
    for ω0 in roots
        ω = solve(LocalDispersionProblem(prob.plasma, k, ω0; closure = prob.closure), alg).omega
        isfinite(ω) || (ω = ω0)
        all(abs(ω - z) > sqrt(alg.atol) for z in polished) && push!(polished, ω)
    end
    return polished
end

function _in_box(region, point = 0)
    ll, ur = region
    return real(ll) ≤ real(point) ≤ real(ur) && imag(ll) ≤ imag(point) ≤ imag(ur)
end

"""
    solve(prob::BranchProblem, alg=ArcLength()) -> DispersionSolution

Track one branch across `prob.ks`. `retcode` is `:Partial` if any `k` failed.
"""
CommonSolve.solve(prob::BranchProblem; kw...) = solve(prob, ArcLength(); kw...)
