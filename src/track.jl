# Arc-length branch continuation
# At each `k` Muller is seeded from a predictor off the previously converged
# roots (linear extrapolation of the last two once both exist, else the last).
# Non-convergence yields NaN for that `k`; the predictor falls back to the last
# finite root so one bad step does not poison the curve. Failed or implausibly
# large jumps trigger a local GRPF search in a square around the predictor.
function _track(plasma, ks, omega0, closure; atol=1e-10, maxiter=100, fallback=true)
    fallback = _track_fallback_config(fallback)
    alg = Muller(; atol, maxiter)
    ωs = ComplexF64[]
    Base.IteratorSize(typeof(ks)) isa Base.HasLength && sizehint!(ωs, length(ks))
    prev = ComplexF64(omega0)   # last finite root (continuation anchor)
    prev2 = ComplexF64(NaN, NaN)
    for k in ks
        # Predictor: linear extrapolation once two finite roots exist.
        guess = isfinite(prev2) && isfinite(prev) ? 2prev - prev2 : prev
        ω = solve(LocalDispersionProblem(plasma, k, guess; closure), alg).omega
        failed = !isfinite(ω)
        jumped = !failed && fallback !== nothing &&
                 _track_jump_too_large(ω, guess, prev, prev2, fallback.jump_factor)
        if fallback !== nothing && (failed || jumped)
            radius = _track_fallback_radius(fallback.radius, guess, prev, prev2,
                                            fallback.radius_factor)
            ωfallback = _track_fallback(plasma, k, guess, radius, fallback.tol, atol, maxiter, closure)
            isfinite(ωfallback) && (ω = ωfallback)
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

const _TRACK_FALLBACK_DEFAULTS = (;
    radius=nothing, tol=1.0e-3, jump_factor=8.0, radius_factor=6.0
)

function _track_fallback_config(fallback)
    fallback === false && return nothing
    return fallback === true ? _TRACK_FALLBACK_DEFAULTS : merge(_TRACK_FALLBACK_DEFAULTS, fallback)
end

function _track_jump_too_large(ω, guess, prev, prev2, factor)
    isfinite(prev) && isfinite(prev2) || return false
    step = abs(prev - prev2)
    step > 0 || return false
    return abs(ω - guess) > factor * step
end

function _track_fallback_radius(radius, guess, prev, prev2, factor)
    isnothing(radius) || return float(radius)
    if isfinite(prev) && isfinite(prev2)
        step = abs(prev - prev2)
        step > 0 && return factor * step
    end
    return max(0.1 * abs(guess), sqrt(eps(Float64)))
end

function _track_fallback(plasma, k, guess, radius, tol, atol, maxiter, closure)
    region = (guess - radius * (1 + im), guess + radius * (1 + im))
    roots = solve(GlobalDispersionProblem(plasma, k, region; closure), GRPF(; tol)).omega
    isempty(roots) && return ComplexF64(NaN, NaN)

    alg = Muller(; atol, maxiter)
    polished = ComplexF64[]
    for root in sort(roots; by=z -> abs(z - guess))
        ω = solve(LocalDispersionProblem(plasma, k, root; closure), alg).omega
        isfinite(ω) || continue
        all(abs(ω - z) > sqrt(atol) for z in polished) && push!(polished, ω)
    end
    isempty(polished) && return roots[argmin(abs.(roots .- guess))]
    return polished[argmin(abs.(polished .- guess))]
end
