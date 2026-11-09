module VlasovMaxwellDispersionRootsExt

using VlasovMaxwellDispersion: _seed_solution, DispersionProblem, DispersionSolution, SolveStats, Seed, Wavevector, residual, _seed_offset, ReturnCode
import VlasovMaxwellDispersion: _solve
import Roots

# find_zero terminates on step size, so it both throws on the iteration cap AND can return a finite step-converged non-root.
# Guard with an explicit residual check so a finite ω genuinely means |f(ω)|≈0
function _solve(prob::DispersionProblem{<:Seed, <:Wavevector}, alg::Roots.AbstractUnivariateZeroMethod; atol=1.0e-10, maxevals=100, kw...)
    t0 = time_ns()
    nevals = Ref(0)
    rawf = prob.f
    f = ω -> (nevals[] += 1; rawf(ω))
    ω0 = prob.target[]
    h = _seed_offset(ω0)
    ω = try
        z = ComplexF64(Roots.find_zero(f, ω0 + h * im, alg; atol, maxevals, kw...))
        abs(f(z)) <= sqrt(atol) * max(abs(f(ω0)), 1) ? z : ComplexF64(NaN, NaN)
    catch
        ComplexF64(NaN, NaN)
    end
    ok = isfinite(ω)
    stats = SolveStats(nevals[], (time_ns() - t0) / 1.0e9)
    return _seed_solution(prob, ω, residual(rawf, ω), stats, ok ? ReturnCode.Success : ReturnCode.Failure, alg)
end

end
