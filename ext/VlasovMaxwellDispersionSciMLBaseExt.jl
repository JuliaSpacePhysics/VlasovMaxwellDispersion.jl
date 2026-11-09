module VlasovMaxwellDispersionSciMLBaseExt

# Lets a point `DispersionProblem` be polished by any SciML initial-guess solver
# (Halley, Newton, Broyden, …) through the same `solve(prob, alg)` verb.
# Caveat: det 𝒟 ∈ ℂ; solvers that assume a real residual
# (Broyden/DFSane compare |f| via `isless`) will error. Halley works.
using VlasovMaxwellDispersion: _seed_solution, DispersionProblem, DispersionSolution, SolveStats, Seed, Wavevector, DispersionFunction,
    residual, prepare, ReturnCode
import VlasovMaxwellDispersion: _solve
import SciMLBase

function _solve(prob::DispersionProblem{<:Seed, <:Wavevector}, alg::SciMLBase.AbstractNonlinearAlgorithm; kwargs...)
    t0 = time_ns()
    nevals = Ref(0)
    f = DispersionFunction(prepare(prob))
    np = SciMLBase.NonlinearProblem((ω, _) -> (nevals[] += 1; f(ω)), complex(prob.target.omega0))
    sol = SciMLBase.solve(np, alg; kwargs...)
    ω = sol.u
    code = if SciMLBase.successful_retcode(sol)
        ReturnCode.Success
    elseif sol.retcode === SciMLBase.ReturnCode.MaxIters
        ReturnCode.MaxIters
    else
        ReturnCode.Failure
    end
    stats = SolveStats(nevals[], (time_ns() - t0) / 1.0e9)
    return _seed_solution(prob, ω, residual(f, ω), stats, code, alg)
end

end
