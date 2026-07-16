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

mutable struct ContinuationCache{P, A, KV, KE, W, R}
    prob::P
    alg::A
    ks::KV
    omega::Vector{W}
    resid::Vector{R}
    seed::W
    khist::Vector{KE}   # accepted points, oldest first, pairwise distinct, ≤ order+1 retained
    whist::Vector{W}
    nevals::Int
end

_saveat(ks) = collect(ks)
_saveat(k::Wavenumber) = [k]
function _saveat(g::ParameterGeometry)
    grids = paramgrids(g)
    length(grids) == 1 || throw(ArgumentError("higher-dimensional sweeps are not supported yet"))
    return map(wavefun(g), only(grids))
end

function CommonSolve.init(prob::DispersionProblem, alg::Continuation)
    p = prepare(prob)
    ks = _saveat(p.k)
    CT = complex(float(typeof(p.omega0)))
    RT = real(CT)
    return ContinuationCache(p, alg, ks, CT[], RT[], CT(p.omega0), similar(ks, 0), CT[], 0)
end

function CommonSolve.step!(cache::ContinuationCache)
    i = length(cache.omega) + 1
    i > length(cache.ks) && return cache
    ω, res = _advance!(cache, cache.ks[i], 0)
    push!(cache.omega, ω)
    push!(cache.resid, res)
    return cache
end

function CommonSolve.solve!(cache::ContinuationCache)
    t0 = time_ns()
    while length(cache.omega) < length(cache.ks)
        step!(cache)
    end
    time = (time_ns() - t0) / 1.0e9
    (; prob, alg, omega, resid) = cache
    return DispersionSolution(
        omega, resid, SolveStats(cache.nevals, time),
        all(isfinite, omega) ? ReturnCode.Success : ReturnCode.Partial, prob, alg
    )
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
    (; prob, alg) = cache
    sol = solve(DispersionProblem(prob.plasma, guess, kt; closure = prob.closure), alg.base)
    cache.nevals += sol.stats.nevals
    return sol.omega, sol.resid
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
# no coordinate has to be carried across the track. Arclength — rather than a
# projection onto the direction of travel — is what keeps the nodes strictly ordered
# even on a path that curves (a swept `theta` is a circular arc) or doubles back.
_dist(a::Wavenumber, b::Wavenumber) = hypot(a.kperp - b.kperp, a.kz - b.kz)
function _arclengths!(s, khist)
    s[1] = zero(eltype(s))
    for i in 2:length(khist)
        s[i] = s[i - 1] + _dist(khist[i - 1], khist[i])
    end
    return s
end

# Newton divided differences through the retained nodes; degree grows with history.
# Scratch is bump-allocated: against a kinetic VDF this is noise, but on a cheap one
# (ColdVDF: a 72 ns det) it would otherwise be the only heap traffic in the loop.
function _predict(cache, kt)
    (; khist, whist) = cache
    isempty(whist) && return cache.seed
    _established(cache) || return last(whist)
    n = length(whist)
    return @no_escape begin
        s = @alloc(typeof(_dist(kt, kt)), n)
        c = @alloc(eltype(whist), n)
        st = _arclengths!(s, khist)[n] + _dist(last(khist), kt)
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

_kmid(a::Wavenumber, b::Wavenumber) = Wavenumber((a.kperp + b.kperp) / 2, (a.kz + b.kz) / 2)
