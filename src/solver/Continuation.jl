"""
    Continuation(; base=Muller(), order=3, reltol=1e-2, abstol=1e-10, maxsubdiv=10)

Track one root branch along an ordered `k` path. Each step extrapolates a degree-`order`
polynomial through the last accepted roots and refines it with `base`.

A step is accepted below `reltol*|ω| + abstol`; otherwise the k-interval is bisected 
and retried, up to `maxsubdiv` levels.

An interval that stays unresolved at `maxsubdiv` yields `NaN` and `ReturnCode.Partial`.
"""
Base.@kwdef struct Continuation{B, T}
    base::B = Muller()
    order::Int = 3
    reltol::T = 1.0e-2
    abstol::T = 1.0e-10
    maxsubdiv::Int = 10
end

mutable struct ContinuationCache{P, A, B, KV, KE, W, R}
    prob::P
    alg::A
    base::B             # `alg.base`'s seeded solver
    ks::KV
    omega::Vector{W}
    resid::Vector{R}
    seed::W
    khist::Vector{KE}   # accepted points, oldest first, pairwise distinct, ≤ order+1 retained
    whist::Vector{W}
    nevals::Int
end

function _saveat(k)
    ks = wavevectors(k)
    ndims(ks) ≤ 1 || throw(ArgumentError("higher-dimensional sweeps are not supported yet"))
    return collect(vec(ks))
end

# Seed anchored mid-path (`at`): track outward both ways and stitch to path order.
function __solve(prob::DispersionProblem{<:Seed}, alg::Continuation; kw...)
    isnothing(prob.target.at) && return _run!(_cache(prob, alg; kw...))
    prob = _prepare(prob, _options(kw))   # once for both directions
    ks = _saveat(prob.k)
    j = argmin(i -> norm(prob.target.at - ks[i]), eachindex(ks))
    fwd = _run!(_cache(_reseed(prob, ks[j:end]), alg; kw...))
    bwd = _run!(_cache(_reseed(prob, reverse(ks[1:j])), alg; kw...))
    return _splice(bwd, fwd, prob, alg)
end
_reseed(prob, ks) = remake(prob; target = Seed(prob.target.omega0), k = ks)

# bwd runs anchor→start over the reversed path; drop its anchor node (fwd owns it).
function _splice(bwd, fwd, prob, alg)
    omega = vcat(reverse(bwd.omega)[1:end-1], fwd.omega)
    resid = vcat(reverse(bwd.resid)[1:end-1], fwd.resid)
    stats = bwd.stats + fwd.stats
    code = successful_retcode(fwd) && successful_retcode(bwd) ? ReturnCode.Success : ReturnCode.Partial
    return DispersionSolution(omega, resid, stats, code, prob, alg)
end

function _cache(prob, alg::Continuation; kw...)
    opts = _options(kw)
    p = _prepare(prob, opts)
    ks = _saveat(p.k)
    ω0 = p.target.omega0
    CT = complex(float(typeof(ω0)))
    RT = real(CT)
    return ContinuationCache(p, alg, seedsolver(p, alg.base; opts...), ks, CT[], RT[], CT(ω0), similar(ks, 0), CT[], 0)
end

function _run!(cache::ContinuationCache)
    t0 = time_ns()
    for k in cache.ks
        ω, res = _advance!(cache, k, 0)
        push!(cache.omega, ω)
        push!(cache.resid, res)
    end
    time = (time_ns() - t0) / 1.0e9
    (; prob, alg, omega, resid) = cache
    code = all(isfinite, omega) ? ReturnCode.Success : ReturnCode.Partial
    return DispersionSolution(omega, resid, SolveStats(cache.nevals, time), code, prob, alg)
end

# Step to `kt`, bisecting the interval from the last accepted point on rejection.
function _advance!(cache, kt, depth)
    guess = _predict(cache, kt)
    ω, res = _refine(cache, kt, guess)
    if _accept(cache, ω, guess)
        _push_accepted!(cache, kt, ω)
        return ω, res
    end
    unresolved = (_complex_nan(guess), oftype(res, NaN))
    stuck = isempty(cache.khist) || vec3(kt) == vec3(last(cache.khist))
    (stuck || depth ≥ cache.alg.maxsubdiv) && return unresolved
    ωm, _ = _advance!(cache, _kmid(last(cache.khist), kt), depth + 1)
    isfinite(ωm) || return unresolved
    return _advance!(cache, kt, depth + 1)
end

function _refine(cache, kt, guess)
    r = cache.base(kt, guess)
    cache.nevals += r.nevals
    return r.omega, r.resid
end

function _push_accepted!(cache, kt, ω)
    # Re-refining the same k replaces its node rather than duplicating it: distinct
    # nodes are what keep the divided-difference denominators non-zero.
    if !isempty(cache.khist) && vec3(kt) == vec3(last(cache.khist))
        cache.whist[end] = ω
        return cache
    end
    push!(cache.khist, kt)
    push!(cache.whist, ω)
    if length(cache.khist) > cache.alg.order + 1
        popfirst!(cache.khist)
        popfirst!(cache.whist)
    end
    return cache
end

# Arclength of the retained nodes along the path travelled, and of `kt` beyond them.
# The origin is arbitrary (here, the oldest retained node): polynomial extrapolation
# is invariant under an affine shift of the abscissa, so only differences matter and
# no coordinate has to be carried across the track.
function _arclengths!(s, khist)
    s[1] = zero(eltype(s))
    for i in 2:length(khist)
        s[i] = s[i - 1] + norm(khist[i] - khist[i - 1])
    end
    return s
end

# Newton divided differences through the retained nodes; degree grows with history.
function _predict(cache, kt)
    (; khist, whist) = cache
    isempty(whist) && return cache.seed
    _established(cache) || return last(whist)
    n = length(whist)
    return @no_escape begin
        s = @alloc(typeof(norm(kt)), n)
        c = @alloc(eltype(whist), n)
        st = _arclengths!(s, khist)[n] + norm(kt - last(khist))
        copyto!(c, whist)
        for j in 2:n, i in n:-1:j
            c[i] = (c[i] - c[i - 1]) / (s[i] - s[i - j + 1])
        end
        v = c[n]
        for i in (n - 1):-1:1
            v = v * (st - s[i]) + c[i]
        end
        v
    end
end

# Two nodes ⇒ an extrapolant exists, and its error gauges the step.
_established(cache) = length(cache.khist) ≥ 2 && isfinite(last(cache.whist))

function _accept(cache, ω, guess)
    isfinite(ω) || return false
    _established(cache) || return true
    scale = max(abs(last(cache.whist)), abs(ω))
    return abs(ω - guess) ≤ cache.alg.reltol * scale + cache.alg.abstol
end

_kmid(a::Wavevector, b::Wavevector) = Wavevector((a.k_perp + b.k_perp) / 2, (a.k_para + b.k_para) / 2)
