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

mutable struct ArcLengthCache{P, A, K, W}
    prob::P
    alg::A
    ks::K
    omega::Vector{W}
    prev::W
    prev2::W
    nevals::Int
end

function CommonSolve.init(prob::DispersionProblem, alg::ArcLength)
    ks = collect(prob.k)
    CT = complex(float(typeof(prob.omega0)))
    prev = CT(prob.omega0)
    return ArcLengthCache(prob, alg, ks, CT[], prev, _complex_nan(prev), 0)
end

function CommonSolve.step!(cache::ArcLengthCache)
    i = length(cache.omega) + 1
    i > length(cache.ks) && return cache
    (; prob, alg, prev, prev2) = cache
    k = cache.ks[i]
    guess = isfinite(prev2) && isfinite(prev) ? 2prev - prev2 : prev
    local_sol = solve(DispersionProblem(prob.plasma, guess, k; closure = prob.closure), alg.base)
    cache.nevals += local_sol.nevals
    ω = local_sol.omega
    if _needs_fallback(alg.fallback, ω, guess, prev, prev2)
        radius = _fallback_radius(alg.fallback, guess, prev, prev2)
        region = (guess - radius * (1 + im), guess + radius * (1 + im))
        survey = solve(GlobalDispersionProblem(prob.plasma, region, k; closure = prob.closure); refine = alg.base)
        cache.nevals += survey.nevals
        roots = [b.omega for b in survey.roots]
        ωfb = isempty(roots) ? complex(NaN) : roots[argmin(abs.(roots .- guess))]
        isfinite(ωfb) && (ω = ωfb)
    end
    push!(cache.omega, ω)
    if isfinite(ω)
        cache.prev2 = prev
        cache.prev = ω
    end
    return cache
end

function CommonSolve.solve!(cache::ArcLengthCache)
    while length(cache.omega) < length(cache.ks)
        step!(cache)
    end
    (; prob, alg, omega) = cache
    res = [residual(prob.plasma, ω, k; closure = prob.closure) for (k, ω) in zip(cache.ks, omega)]
    return DispersionSolution(
        omega, res, cache.nevals,
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
