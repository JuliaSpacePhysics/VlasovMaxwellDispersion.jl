module VlasovMaxwellDispersionRootsExt

using VlasovMaxwellDispersion: LocalDispersionProblem, DispersionSolution, residual
import CommonSolve: solve
import Roots

# find_zero terminates on step size, so it both throws on the iteration cap AND can return a finite step-converged non-root.
# Guard with an explicit residual check so a finite ω genuinely means |f(ω)|≈0 (track.jl's fallback relies on)
function solve(prob::LocalDispersionProblem, alg::Roots.AbstractUnivariateZeroMethod; atol=1.0e-10, maxevals=100, kw...)
    f = prob.f
    h = 1.0e-3 * max(abs(prob.omega0), 1.0)
    ω = try
        z = ComplexF64(Roots.find_zero(f, prob.omega0 + h * im, alg; atol, maxevals, kw...))
        abs(f(z)) <= sqrt(atol) * max(abs(f(prob.omega0)), 1) ? z : ComplexF64(NaN, NaN)
    catch
        ComplexF64(NaN, NaN)
    end
    ok = isfinite(ω)
    return DispersionSolution(ω, nothing, residual(prob, ω), ok ? :Success : :Failure, prob, alg)
end

end
