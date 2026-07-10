"""
    Muller(; atol=1e-10, maxiter=100)

Derivative-free complex root polish.
"""
Base.@kwdef struct Muller
    atol::Float64 = 1.0e-10
    maxiter::Int = 100
end

# Absolute-floored relative step capped at 0.1|ω0|.
# ω0===0 has no scale ⇒ small absolute step.
_seed_offset(ω0) = iszero(ω0) ? 1.0e-3 : min(1.0e-3 * max(abs(ω0), 1.0), 0.1 * abs(ω0))

function CommonSolve.solve(prob::DispersionProblem{<:Any, <:Wavenumber}, alg::Muller)
    f = prob.f
    h = _seed_offset(prob.omega0)
    ω = muller(f, prob.omega0 - h, prob.omega0 + h, prob.omega0 + h * im; alg.atol, alg.maxiter)
    ok = isfinite(ω)
    return DispersionSolution(ω, residual(prob, ω), ok ? :Success : :Failure, prob, alg)
end


"""
    muller(f, x0, x1, x2; atol=1e-10, maxiter=100)

Complex-valued Muller's method. Fits a quadratic through `(x0,f0),(x1,f1),(x2,f2)`
and advances to the root closer to `x2`, picking the denominator branch that
maximizes `|denom|` for numerical stability. 
    
Returns `NaN+NaN*im` if `maxiter` is exhausted without reaching `atol` on the step size or `|f|`.
"""
function muller(f, x0, x1, x2; atol = 1.0e-10, maxiter = 100)
    x0, x1, x2 = promote(complex(float(x0)), complex(float(x1)), complex(float(x2)))
    nan = typeof(x2)(NaN, NaN)
    f0, f1, f2 = f(x0), f(x1), f(x2)
    for _ in 1:maxiter
        h1 = x1 - x0
        h2 = x2 - x1
        if h1 == 0 || h2 == 0
            return nan
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
            return nan
        end

        dx = -2c / denom
        x3 = x2 + dx
        f3 = f(x3)
        # Contract trials that hit the determinant's overflow back toward finite f(x2)
        flimit = 1.0e6 * max(abs(f0), abs(f1), abs(f2))
        for _ in 1:12
            (isfinite(f3) && abs(f3) <= flimit) && break
            dx /= 2
            x3 = x2 + dx
            f3 = f(x3)
        end
        isfinite(f3) || return nan # whole neighborhood is bad — bail
        fscale = max(abs(f0), abs(f1), abs(f2), 1)

        if abs(f3) <= atol * fscale ||
                (abs(dx) <= atol * max(abs(x3), 1) && abs(f3) <= sqrt(atol) * fscale)
            return x3
        end

        x0, x1, x2 = x1, x2, x3
        f0, f1, f2 = f1, f2, f3
    end
    return abs(f2) <= atol ? x2 : nan
end
