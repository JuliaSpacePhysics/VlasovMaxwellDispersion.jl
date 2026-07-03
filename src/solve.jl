dispersion_function(prob::AbstractDispersionProblem) =
    ω -> det(dispersion_tensor(prob.plasma, ω, prob.k; closure = prob.closure))

@inline Base.getproperty(prob::LocalDispersionProblem, s::Symbol) =
    s === :f ? _scaled_dispersion_function(prob) : getfield(prob, s)

# GRPF is phase-based and has no seed. det𝒟 has a genuine pole at ω=0 (Maxwell
# light term c²k²/ω², order ≤6) whose winding partially cancels nearby roots, so
# the argument principle sees det(ω²𝒟) = ω⁶det𝒟, which is entire there. The
# deflation leaves an artifact zero of order 6−p at exactly ω=0; solve filters it.
@inline Base.getproperty(prob::GlobalDispersionProblem, s::Symbol) =
    s === :f ? _deflated_dispersion_function(prob) : getfield(prob, s)

function _deflated_dispersion_function(prob::GlobalDispersionProblem)
    f = dispersion_function(prob)
    return ω -> ω^6 * f(ω)
end

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
    solve(prob::GlobalDispersionProblem, alg=GRPF()) -> DispersionSolution

All roots and poles in `prob.region` via the argument principle, applied to the
deflated determinant `det(ω²𝒟) = ω⁶ det𝒟` so the ω=0 light-term pole cannot
cancel or displace nearby roots. Roots within `alg.tol` of ω=0 are dropped:
at that distance they are indistinguishable from the deflation artifact.
"""
function CommonSolve.solve(prob::GlobalDispersionProblem, alg = GRPF())
    roots, poles = _grpf_roots(prob.f, prob.region; alg.tol, alg.params)
    filter!(ω -> abs(ω) > alg.tol, roots)   # deflation artifact at ω=0
    # GRPF roots are mesh-accurate (~tol), so these sit well above ε
    res = [residual(prob, ω) for ω in roots]
    return DispersionSolution(roots, poles, res, isempty(roots) ? :Failure : :Success, prob, alg)
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
