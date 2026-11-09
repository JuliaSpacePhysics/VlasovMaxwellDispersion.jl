module VlasovMaxwellDispersionSciMLBaseExt

# Lets a point `DispersionProblem` be polished by any SciML initial-guess solver
# (Halley, Newton, Broyden, …) through the same `solve(prob, alg)` verb.
# Caveat: det 𝒟 ∈ ℂ; solvers that assume a real residual
# (Broyden/DFSane compare |f| via `isless`) will error. Halley works.
using VlasovMaxwellDispersion: _seedsolver, _split_options, ReturnCode
import VlasovMaxwellDispersion: findroot
import PlasmaDispersionBase: seedsolver
import SciMLBase

seedsolver(prob, alg::SciMLBase.AbstractNonlinearAlgorithm; kw...) = _seedsolver(prob, alg, _split_options(kw)...)

function findroot(alg::SciMLBase.AbstractNonlinearAlgorithm, D, ω0; kw...)
    nevals = Ref(0)
    np = SciMLBase.NonlinearProblem((ω, _) -> (nevals[] += 1; D(ω)), complex(ω0))
    sol = SciMLBase.solve(np, alg; kw...)
    code = SciMLBase.successful_retcode(sol) ? ReturnCode.Success :
        sol.retcode === SciMLBase.ReturnCode.MaxIters ? ReturnCode.MaxIters : ReturnCode.Failure
    return sol.u, nevals[], code
end

end
