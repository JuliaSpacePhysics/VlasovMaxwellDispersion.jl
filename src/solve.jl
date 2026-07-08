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
residual(prob, ω) =
    residual(prob.plasma, ω, prob.k; closure = prob.closure)

function _scaled_dispersion_function(prob::LocalDispersionProblem)
    f = dispersion_function(prob)
    s = _hadamard(dispersion_tensor(prob.plasma, prob.omega0, prob.k; closure = prob.closure))
    return isfinite(s) && s > 0 ? (ω -> f(ω) / s) : f
end

include("solver/muller.jl")

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
    solve(prob::GlobalDispersionProblem, alg=GRPF()) -> DispersionSolution

All roots and poles in `prob.region` via the argument principle.
"""
function CommonSolve.solve(prob::GlobalDispersionProblem, alg = GRPF())
    roots, poles = _grpf_roots(prob.f, prob.region; alg.tol, alg.params)
    # Drop roots the mesh cannot separate,
    # The artifact sits at |ω| ≲ tol (mesh accuracy); 2·tol gives a one-cell margin.
    _in_box(prob.region) && filter!(ω -> abs(ω) > 2alg.tol, roots)
    res = [residual(prob, ω) for ω in roots]
    return DispersionSolution(roots, poles, res, isempty(roots) ? :Failure : :Success, prob, alg)
end

function _in_box(region, point = 0)
    ll, ur = region
    return real(ll) ≤ real(point) ≤ real(ur) && imag(ll) ≤ imag(point) ≤ imag(ur)
end

# Low-level GRPF over a complex box; returns (roots, poles).
function _grpf_roots(f, region; tol = 1.0e-3, params = nothing)
    lowerleft, upperright = ComplexF64(region[1]), ComplexF64(region[2])
    origcoords = rectangulardomain(lowerleft, upperright, tol)
    p = @something params GRPFParams(5000, tol, false)
    roots, poles = grpf(f, origcoords, p)
    return ComplexF64.(roots), ComplexF64.(poles)
end

"""
    solve(prob::BranchProblem, alg=ArcLength()) -> DispersionSolution

Track one branch across `prob.ks`. `retcode` is `:Partial` if any `k` failed.
"""
function CommonSolve.solve(prob::BranchProblem, alg = ArcLength())
    ωs = _track(prob.plasma, prob.ks, prob.omega0, prob.closure; alg.atol, alg.maxiter, alg.fallback)
    res = [residual(prob.plasma, ω, k; closure = prob.closure) for (k, ω) in zip(prob.ks, ωs)]
    return DispersionSolution(ωs, nothing, res, all(isfinite, ωs) ? :Success : :Partial, prob, alg)
end
