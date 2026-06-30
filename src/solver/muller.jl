"""
    Muller(; atol=1e-10, maxiter=100)

Derivative-free complex root polish for [`LocalDispersionProblem`](@ref).
"""
Base.@kwdef struct Muller <: DispersionAlgorithm
    atol::Float64 = 1.0e-10
    maxiter::Int = 100
end

# Three seeds are clustered around `omega0` (small relative+absolute perturbation).
function CommonSolve.solve(prob::LocalDispersionProblem, alg::Muller)
    f = residual(prob)
    h = 1.0e-3 * max(abs(prob.omega0), 1.0)
    ω = muller(f, prob.omega0 - h, prob.omega0 + h, prob.omega0 + h * im; alg.atol, alg.maxiter)
    ok = isfinite(ω)
    return DispersionSolution(ω, nothing, ok ? abs(f(ω)) : NaN, ok ? :Success : :Failure, prob, alg)
end


"""
    muller(f, x0, x1, x2; atol=1e-10, maxiter=100)

Complex-valued Muller's method. Fits a quadratic through `(x0,f0),(x1,f1),(x2,f2)`
and advances to the root closer to `x2`, picking the denominator branch that
maximizes `|denom|` for numerical stability. 
    
Returns `NaN+NaN*im` if `maxiter` is exhausted without reaching `atol` on the step size or `|f|`.
"""
function muller(f, x0, x1, x2; atol=1.0e-10, maxiter=100)
    x0, x1, x2 = ComplexF64(x0), ComplexF64(x1), ComplexF64(x2)
    f0, f1, f2 = f(x0), f(x1), f(x2)
    # The dispersion determinant overflows to Inf/NaN for strongly damped ω (the VDF's analytic
    # continuation is evaluated far off the real velocity axis). Without a guard the quadratic
    # step then produces NaN and the iteration wanders into — and stalls on — that region. Bail
    # if the seeds are already non-finite; below, a non-finite trial step is contracted back
    # toward the last finite point instead of being accepted.
    (isfinite(f0) && isfinite(f1) && isfinite(f2)) || return ComplexF64(NaN, NaN)
    for _ in 1:maxiter
        h1 = x1 - x0
        h2 = x2 - x1
        if h1 == 0 || h2 == 0
            return ComplexF64(NaN, NaN)
        end
        delta1 = (f1 - f0) / h1
        delta2 = (f2 - f1) / h2
        a = (delta2 - delta1) / (h2 + h1)
        b = a * h2 + delta2
        c = f2

        disc = sqrt(b^2 - 4a * c)
        denom_plus = b + disc
        denom_minus = b - disc
        denom = abs(denom_plus) > abs(denom_minus) ? denom_plus : denom_minus
        if denom == 0
            return ComplexF64(NaN, NaN)
        end

        dx = -2c / denom
        x3 = x2 + dx
        f3 = f(x3)
        # Contract a non-finite / wildly-diverging trial back toward x2 (which is finite).
        nc = 0
        while (!isfinite(f3) || abs(f3) > 1.0e6 * max(abs(f2), 1)) && nc < 12
            dx /= 2
            x3 = x2 + dx
            f3 = f(x3)
            nc += 1
        end
        isfinite(f3) || return ComplexF64(NaN, NaN)
        fscale = max(abs(f0), abs(f1), abs(f2), 1)

        if abs(f3) <= atol * fscale ||
           (abs(dx) <= atol * max(abs(x3), 1) && abs(f3) <= sqrt(atol) * fscale)
            return x3
        end

        x0, x1, x2 = x1, x2, x3
        f0, f1, f2 = f1, f2, f3
    end
    return abs(f2) <= atol ? x2 : ComplexF64(NaN, NaN)
end