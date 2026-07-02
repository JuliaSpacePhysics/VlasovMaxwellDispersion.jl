module VlasovMaxwellDispersionSciMLBaseExt

# Lets a `LocalDispersionProblem` be polished by any SciML initial-guess solver
# (Halley, Newton, Broyden, …) through the same `solve(prob, alg)` verb.
# Caveat: det 𝒟 ∈ ℂ; solvers that assume a real residual 
# (Broyden/DFSane compare |f| via `isless`) will error. Halley works.
using VlasovMaxwellDispersion: LocalDispersionProblem, DispersionSolution, residual
import CommonSolve: solve
import SciMLBase

function solve(prob::LocalDispersionProblem, alg::SciMLBase.AbstractNonlinearAlgorithm; kwargs...)
    np = SciMLBase.NonlinearProblem((ω, _) -> prob.f(ω), prob.omega0)
    sol = solve(np, alg; kwargs...)
    ω = sol.u
    ok = SciMLBase.successful_retcode(sol)
    return DispersionSolution(ω, nothing, residual(prob, ω), ok ? :Success : :Failure, prob, alg)
end

end
