residual(prob::Union{LocalDispersionProblem,GlobalDispersionProblem}) =
    ω -> det(dispersion_tensor(prob.plasma, ω, prob.k; closure=prob.closure))

"""
    dispersion_residual(plasma, ω, k::Wavenumber; closure=HarmonicSum()) -> Float64

Scale-free "rootness" of `𝒟(ω,k)`: `σ_min(𝒟)/σ_max(𝒟)` ∈ `[0,1]`, zero at a mode.

Prefer this over `abs(det(𝒟))` for convergence, root validation, and cross-solver checks.
`det(𝒟)` floors at `‖𝒟‖³·ε` from catastrophic cancellation — e.g. `|det|` bottoms out around
`1` even at a clean root when the tensor entries are `~1e5` — so an absolute `|det| < tol` test
is scale-dependent and may never trigger. The singular-value ratio is normalized and reaches
`~machine-ε` at a genuine root regardless of the tensor's overall scale.
"""
function dispersion_residual(plasma, ω, k::Wavenumber; kwargs...)
    s = svdvals(Matrix(dispersion_tensor(plasma, ω, k; kwargs...)))
    return s[end] / s[1]
end

include("solver/secant.jl")
include("solver/muller.jl")

"""
    solve(prob::LocalDispersionProblem, alg) -> DispersionSolution

Default alg = Muller().
"""
CommonSolve.solve(prob::LocalDispersionProblem) =
    CommonSolve.solve(prob, Muller())

"""
    solve(prob::GlobalDispersionProblem, alg=GRPF()) -> DispersionSolution

All roots and poles in `prob.region` via the argument principle.
"""
function CommonSolve.solve(prob::GlobalDispersionProblem, alg=GRPF())
    f = residual(prob)
    roots, poles = _grpf_roots(f, prob.region; alg.tol, alg.params)
    return DispersionSolution(roots, poles, nothing, isempty(roots) ? :Failure : :Success, prob, alg)
end

# Low-level GRPF over a complex box; returns (roots, poles).
function _grpf_roots(f, region; tol=1.0e-3, params=nothing)
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
function CommonSolve.solve(prob::BranchProblem, alg=ArcLength())
    ωs = _track(prob.plasma, prob.ks, prob.omega0, prob.closure; alg.atol, alg.maxiter, alg.fallback)
    return DispersionSolution(ωs, nothing, nothing, all(isfinite, ωs) ? :Success : :Partial, prob, alg)
end
