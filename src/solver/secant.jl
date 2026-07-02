"""
    Secant(; atol=1e-10, maxiter=100)

Secant method via `Roots.Order1`. 

Lower convergence order than [`Muller`](@ref) (≈1.62 vs 1.84)
so typically a touch slower; seeded with a small imaginary nudge to leave the real axis.
"""
Base.@kwdef struct Secant <: DispersionAlgorithm
    atol::Float64 = 1.0e-10
    maxiter::Int = 100
end

# find_zero terminates on step size, so it both throws on the iteration cap AND can return a finite step-converged non-root.
# Guard with an explicit residual check so a finite ω genuinely means |f(ω)|≈0 (track.jl's fallback relies on)
function CommonSolve.solve(prob::LocalDispersionProblem, alg::Secant)
    f = prob.f
    h = 1.0e-3 * max(abs(prob.omega0), 1.0)
    ω = try
        z = ComplexF64(find_zero(f, prob.omega0 + h * im, Roots.Order1(); atol=alg.atol, maxevals=alg.maxiter))
        abs(f(z)) <= sqrt(alg.atol) * max(abs(f(prob.omega0)), 1) ? z : ComplexF64(NaN, NaN)
    catch
        ComplexF64(NaN, NaN)
    end
    ok = isfinite(ω)
    return DispersionSolution(ω, nothing, residual(prob, ω), ok ? :Success : :Failure, prob, alg)
end