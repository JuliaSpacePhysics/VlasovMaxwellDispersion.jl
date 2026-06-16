residual(prob::Union{LocalDispersionProblem,GlobalDispersionProblem}) =
    ω -> det(dispersion_tensor(prob.plasma, ω, prob.k; closure=prob.closure))

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
