"""
    JumpFallback(; radius=nothing, jump_factor=8.0, radius_factor=6.0)

Local recovery for [`ArcLength`](@ref), triggered when a step fails to converge
or lands farther than `jump_factor`× the previous step from the predictor: a
global survey in a box of half-width `radius` (default `radius_factor`× the last
step) around the predictor. Keeps the one nearest the predictor. A step that
cannot be recovered becomes `NaN`.
"""
Base.@kwdef struct JumpFallback
    radius::Union{Nothing, Float64} = nothing
    jump_factor::Float64 = 8.0
    radius_factor::Float64 = 6.0
end

"""
    ArcLength(; base=Muller(), fallback=JumpFallback())

Continuing `omega0` by extrapolating from previous converged roots and refining with `base`.
`fallback` ([`JumpFallback`](@ref), or `nothing` to disable) catches failed or
implausibly large jumps; an unrecoverable step becomes `NaN` (`:Partial`) rather
than a silent branch switch.
"""
Base.@kwdef struct ArcLength{B, F}
    base::B = Muller()
    fallback::F = JumpFallback()
end

function CommonSolve.solve(prob::BranchProblem, alg::ArcLength)
    ωs = _track(prob.plasma, prob.ks, prob.omega0, prob.closure, alg)
    res = [residual(prob.plasma, ω, k; closure = prob.closure) for (k, ω) in zip(prob.ks, ωs)]
    return DispersionSolution(ωs, res, all(isfinite, ωs) ? :Success : :Partial, prob, alg)
end

# Predictor falls back to the last finite root so one bad step does not poison the curve.
# Failed or implausibly large jumps trigger `alg.fallback`, a search around the predictor.
function _track(plasma, ks, omega0, closure, alg)
    CT = complex(float(typeof(omega0)))
    ωs = CT[]
    Base.IteratorSize(typeof(ks)) isa Base.HasLength && sizehint!(ωs, length(ks))
    prev = CT(omega0)   # last finite root (continuation anchor)
    prev2 = _complex_nan(prev)
    for k in ks
        # Predictor: linear extrapolation once two finite roots exist.
        guess = isfinite(prev2) && isfinite(prev) ? 2prev - prev2 : prev
        ω = solve(LocalDispersionProblem(plasma, k, guess; closure), alg.base).omega
        if _needs_fallback(alg.fallback, ω, guess, prev, prev2)
            radius = _fallback_radius(alg.fallback, guess, prev, prev2)
            region = (guess - radius * (1 + im), guess + radius * (1 + im))
            survey = solve(GlobalDispersionProblem(plasma, region, k; closure); refine = alg.base)
            roots = [b.omega for b in survey.roots]
            ωfb = isempty(roots) ? complex(NaN) : roots[argmin(abs.(roots .- guess))]
            isfinite(ωfb) && (ω = ωfb)
        end
        push!(ωs, ω)
        if isfinite(ω)
            prev2 = prev
            prev = ω
        end
        # Diverged step keeps prev/prev2 (last good anchor) so we retry from it.
    end
    return ωs
end

_needs_fallback(::Nothing, ω, guess, prev, prev2) = false
function _needs_fallback(fb, ω, guess, prev, prev2)
    isfinite(ω) || return true
    isfinite(prev) && isfinite(prev2) || return false
    step = abs(prev - prev2)
    return step > 0 && abs(ω - guess) > fb.jump_factor * step
end

function _fallback_radius(fb, guess, prev, prev2)
    isnothing(fb.radius) || return fb.radius
    if isfinite(prev) && isfinite(prev2)
        step = abs(prev - prev2)
        step > 0 && return fb.radius_factor * step
    end
    return max(0.1 * abs(guess), sqrt(eps(typeof(abs(guess)))))
end
