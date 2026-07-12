# Generic global survey: per-slice zero finding + cross-slice linking. Any alg
# implementing `discover(alg, f, region) -> (zeros, saturated)` and
# `_origin_gate(alg, diag)` solves GlobalDispersionProblem at any swept
# dimension through this one path.

function CommonSolve.init(
        prob::GlobalDispersionProblem, alg;
        refine = Muller(), linking = nothing
    )
    ll, ur = prob.region
    diag = hypot(real(ur - ll), imag(ur - ll))
    linking = @something linking (; maxgap = 1, gate = diag / 8)

    return SurveyCache(prepare(prob), alg, refine, linking)
end

polish!(f, ωs, ::Nothing) = ωs, 0

"""
    link(alg, grids, vals)

Link per-grid-point complex value sets `vals[I]` into continuous sheets over
the m-D parameter grid `grids`. 
Values match extrapolated sheet tips greedily within `gate`, bridging up to
`maxgap` empty points; unmatched values start new sheets. Extrapolation (not
plain nearest-tip) is what keeps two crossing branches from swapping partners.
"""
function link(alg, grids, vals)
    m = length(grids)
    P = float(promote_type(map(eltype, grids)...))
    V = ComplexF64
    sheets = Vector{Tuple{NTuple{m, P}, V}}[]
    tips = Dict{CartesianIndex{m}, Vector{Tuple{Int, V, V}}}()  # point ⇒ (sheet id, value, velocity)
    for c in CartesianIndices(size(vals))
        p = ntuple(j -> P(grids[j][c[j]]), m)
        vs = V.(vals[c])
        ref = Dict{Int, Tuple{V, V, Int}}()          # sid ⇒ (tip, velocity, gap l); closest l wins
        for j in 1:m, l in 1:(alg.maxgap + 1)
            c[j] - l ≥ 1 || break
            nb = c - CartesianIndex(ntuple(i -> i == j ? l : 0, m))
            for (sid, v, dv) in get(tips, nb, Tuple{Int, V, V}[])
                (haskey(ref, sid) && ref[sid][3] ≤ l) || (ref[sid] = (v, dv, l))
            end
        end
        cand = sort!(
            [
                (abs(vs[i] - (v + l * dv)), i, sid, v, l)
                    for (sid, (v, dv, l)) in ref for i in eachindex(vs)
            ]; by = first
        )
        usedv = falses(length(vs))
        useds = Set{Int}()
        assigned = Tuple{Int, V, V}[]
        for (d, i, sid, v, l) in cand
            (d ≤ alg.gate && !usedv[i] && sid ∉ useds) || continue
            usedv[i] = true
            push!(useds, sid)
            push!(sheets[sid], (p, vs[i]))
            push!(assigned, (sid, vs[i], (vs[i] - v) / l))
        end
        for i in eachindex(vs)
            usedv[i] && continue
            push!(sheets, [(p, vs[i])])
            push!(assigned, (length(sheets), vs[i], zero(V)))
        end
        tips[c] = assigned
    end
    return sheets
end

struct SurveyCache{P, A, R, L}
    prob::P
    alg::A
    refine::R
    linking::L
end

function _tforeach(f, xs)
    Threads.@threads for x in xs
        f(x)
    end
    return nothing
end

function CommonSolve.solve!(cache::SurveyCache)
    t0 = time_ns()
    (; prob, alg) = cache
    grids = paramgrids(prob.geometry)
    m = length(grids)
    kf = wavefun(prob.geometry)
    ll, ur = prob.region
    diag = hypot(real(ur - ll), imag(ur - ll))
    gate0 = _in_box(prob.region) ? _origin_gate(alg, diag) : 0.0
    dims = map(length, grids)
    zv = Array{Vector{ComplexF64}}(undef, dims)
    nev = zeros(Int, dims)
    _tforeach(CartesianIndices(dims)) do c
        p = ntuple(j -> grids[j][c[j]], m)
        k = kf(p...)
        f = ω -> (det(wave_dispersion_tensor(prob.plasma, ω, k; closure = prob.closure)))
        zs, nevals = discover(alg, f, (ll, ur))
        polished, np = polish!(f, zs, cache.refine)
        checked = filter!(polished) do z
            isfinite(z) && _in_box((ll, ur), z) && abs(z) > gate0
        end
        zs = consolidate(checked; atol = 1.0e-4 * diag)
        zv[c], nev[c] = zs, nevals + np
        return nothing
    end
    stats = SolveStats(sum(nev), (time_ns() - t0) / 1.0e9)
    return build_solution(cache, grids, zv, stats)
end

function build_solution(cache::SurveyCache, grids, values, stats)
    (; prob, alg) = cache
    m = length(grids)
    kf = wavefun(prob.geometry)
    if m == 0
        k0 = kf()
        roots = [
            DispersionBranch(ω, k0, residual(prob, ω, k0))
                for ω in values[] if _in_box(prob.region, ω)
        ]
        return SurveySolution(roots, stats, _retcode(roots), prob, alg)
    end
    sheets = link(cache.linking, grids, values)
    filter!(sh -> any(_in_box(prob.region, ω) for (_, ω) in sh), sheets)
    T = _realtype(prob)
    makebranch(sh) = begin
        ωs = [Complex{T}(ω) for (_, ω) in sh]
        ks = [kf(p...) for (p, _) in sh]
        res = [residual(prob.plasma, ω, k; closure = prob.closure) for (ω, k) in zip(ωs, ks)]
        DispersionBranch(ωs, ks, res)
    end
    roots = map(makebranch, sheets)
    return SurveySolution(roots, stats, _retcode(roots), prob, alg)
end

_retcode(roots) = isempty(roots) ? :Failure : :Success

# Deduplicate candidate roots using absolute distance `atol`
function consolidate(points; atol)
    out = empty(points)
    for point in points
        any(other -> abs(point - other) < atol, out) || push!(out, point)
    end
    return out
end
