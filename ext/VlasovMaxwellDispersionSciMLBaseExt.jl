module VlasovMaxwellDispersionSciMLBaseExt

# Lets a point `DispersionProblem` be polished by any SciML initial-guess solver
# (Halley, Newton, Broyden, …) through the same `solve(prob, alg)` verb.
# Caveat: det 𝒟 ∈ ℂ; solvers that assume a real residual 
# (Broyden/DFSane compare |f| via `isless`) will error. Halley works.
using VlasovMaxwellDispersion: DispersionProblem, DispersionSolution, SolveStats, Wavenumber, residual
import CommonSolve: solve
import SciMLBase

function solve(prob::DispersionProblem{<:Any, <:Wavenumber}, alg::SciMLBase.AbstractNonlinearAlgorithm; kwargs...)
    t0 = time_ns()
    nevals = Ref(0)
    np = SciMLBase.NonlinearProblem((ω, _) -> (nevals[] += 1; prob.f(ω)), complex(prob.omega0))
    sol = solve(np, alg; kwargs...)
    ω = sol.u
    ok = SciMLBase.successful_retcode(sol)
    stats = SolveStats(nevals[], (time_ns() - t0) / 1.0e9)
    return DispersionSolution(ω, residual(prob, ω), stats, ok ? :Success : :Failure, prob, alg)
end

end
