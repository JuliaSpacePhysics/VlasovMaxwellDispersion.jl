module VlasovMaxwellDispersionSciMLBaseExt

# Lets a point `DispersionProblem` be polished by any SciML initial-guess solver
# (Halley, Newton, Broyden, …) through the same `solve(prob, alg)` verb.
# Caveat: det 𝒟 ∈ ℂ; solvers that assume a real residual 
# (Broyden/DFSane compare |f| via `isless`) will error. Halley works.
using VlasovMaxwellDispersion: DispersionProblem, DispersionSolution, Wavenumber, residual
import CommonSolve: solve
import SciMLBase

function solve(prob::DispersionProblem{<:Any, <:Wavenumber}, alg::SciMLBase.AbstractNonlinearAlgorithm; kwargs...)
    np = SciMLBase.NonlinearProblem((ω, _) -> prob.f(ω), prob.omega0)
    sol = solve(np, alg; kwargs...)
    ω = sol.u
    ok = SciMLBase.successful_retcode(sol)
    return DispersionSolution(ω, residual(prob, ω), ok ? :Success : :Failure, prob, alg)
end

end
