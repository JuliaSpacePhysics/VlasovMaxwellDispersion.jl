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

mutable struct ArcLengthCache{P, A, K, W, R}
    prob::P
    alg::A
    ks::K
    omega::Vector{W}
    res::Vector{R}
    prev::W
    prev2::W
    nevals::Int
end

function CommonSolve.init(prob::DispersionProblem, alg::ArcLength)
    p = prepare(prob)
    CT = complex(float(typeof(p.omega0)))
    prev = CT(p.omega0)
    return ArcLengthCache(p, alg, collect(p.k), CT[], real(CT)[], prev, _complex_nan(prev), 0)
end

function CommonSolve.step!(cache::ArcLengthCache)
    i = length(cache.omega) + 1
    i > length(cache.ks) && return cache
    (; prob, alg, prev, prev2) = cache
    k = cache.ks[i]
    guess = isfinite(prev2) && isfinite(prev) ? 2prev - prev2 : prev
    local_sol = solve(DispersionProblem(prob.plasma, guess, k; closure = prob.closure), alg.base)
    cache.nevals += local_sol.stats.nevals
    ω, res = local_sol.omega, local_sol.resid
    if _needs_fallback(alg.fallback, ω, guess, prev, prev2)
        radius = _fallback_radius(alg.fallback, guess, prev, prev2)
        region = (guess - radius * (1 + im), guess + radius * (1 + im))
        survey = solve(GlobalDispersionProblem(prob.plasma, region, k; closure = prob.closure); refine = alg.base)
        cache.nevals += survey.stats.nevals
        if !isempty(survey.roots)
            b = argmin(b -> abs(b.omega - guess), survey.roots)
            isfinite(b.omega) && ((ω, res) = (b.omega, b.resid))
        end
    end
    push!(cache.omega, ω)
    push!(cache.res, res)
    if isfinite(ω)
        cache.prev2 = prev
        cache.prev = ω
    end
    return cache
end

function CommonSolve.solve!(cache::ArcLengthCache)
    t0 = time_ns()
    while length(cache.omega) < length(cache.ks)
        step!(cache)
    end
    time = (time_ns() - t0) / 1.0e9
    (; prob, alg, omega, res) = cache
    return DispersionSolution(
        omega, res, SolveStats(cache.nevals, time),
        all(isfinite, omega) ? :Success : :Partial, prob, alg
    )
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
