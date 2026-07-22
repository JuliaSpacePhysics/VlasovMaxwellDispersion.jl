module VlasovMaxwellDispersionRootsExt

using VlasovMaxwellDispersion: DispersionProblem, DispersionSolution, SolveStats, Seed, Wavenumber, residual, _seed_offset, ReturnCode
import CommonSolve: solve
import Roots

# find_zero terminates on step size, so it both throws on the iteration cap AND can return a finite step-converged non-root.
# Guard with an explicit residual check so a finite ω genuinely means |f(ω)|≈0
function solve(prob::DispersionProblem{<:Seed, <:Wavenumber}, alg::Roots.AbstractUnivariateZeroMethod; atol=1.0e-10, maxevals=100, kw...)
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
    return DispersionSolution(ω, residual(prob, ω), stats, ok ? ReturnCode.Success : ReturnCode.Failure, prob, alg)
end

end
