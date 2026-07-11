# Generic global survey: per-slice zero finding + cross-slice linking. Any alg
# implementing `_slice_zeros(alg, f, region) -> (zeros, saturated)` and
# `_origin_gate(alg, diag)` solves GlobalDispersionProblem at any swept
# dimension through this one path.

"""
    link_sheets(grids, vals; gate, maxgap=1) -> Vector{Vector{Tuple{NTuple,ComplexF64}}}

Link per-grid-point complex value sets `vals[I]` into continuous sheets over
the m-D parameter grid `grids`. 
Values match extrapolated sheet tips greedily within `gate`, bridging up to
`maxgap` empty points; unmatched values start new sheets. Extrapolation (not
plain nearest-tip) is what keeps two crossing branches from swapping partners.
"""
function link_sheets(grids, vals; gate, maxgap = 0)
    m = length(grids)
    P = float(promote_type(map(eltype, grids)...))
    V = ComplexF64
    sheets = Vector{Tuple{NTuple{m, P}, V}}[]
    tips = Dict{CartesianIndex{m}, Vector{Tuple{Int, V, V}}}()  # point ⇒ (sheet id, value, velocity)
    for c in CartesianIndices(size(vals))
        p = ntuple(j -> P(grids[j][c[j]]), m)
        vs = V.(vals[c])
        ref = Dict{Int, Tuple{V, V, Int}}()          # sid ⇒ (tip, velocity, gap l); closest l wins
        for j in 1:m, l in 1:(maxgap + 1)
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
            (d ≤ gate && !usedv[i] && sid ∉ useds) || continue
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

struct SurveyCache{P, A, R}
    prob::P
    alg::A
    refine::R
    pad::Float64
    gate::Union{Nothing, Float64}
    maxgap::Int
end

CommonSolve.init(
    prob::GlobalDispersionProblem, alg;
    refine = Muller(), pad = 0.15, gate = nothing, maxgap = 1
) = SurveyCache(prob, alg, refine, pad, gate, maxgap)

_tmap(f, xs) = map(fetch, [Threads.@spawn(f(x)) for x in xs])

function CommonSolve.solve!(cache::SurveyCache)
    t0 = time_ns()
    (; prob, alg, refine, pad, maxgap) = cache
    grids = paramgrids(prob.geometry)
    m = length(grids)
    kf = wavefun(prob.geometry)
    ll, ur = prob.region
    # Soft window: track `pad` past every edge so an edge-dipping branch stays
    # one sheet; only in-box results are returned.
    pll, pur = ll - pad * (ur - ll), ur + pad * (ur - ll)
    diag = hypot(real(ur - ll), imag(ur - ll))
    gate0 = _in_box(prob.region) ? _origin_gate(alg, diag) : 0.0
    dims = map(length, grids)
    zv = Array{Vector{ComplexF64}}(undef, dims)
    nev = zeros(Int, dims)
    sat = fill(false, dims)
    _tmap(CartesianIndices(dims)) do c
        p = ntuple(j -> grids[j][c[j]], m)
        k = kf(p...)
        n = Ref(0)
        f = ω -> (n[] += 1; det(wave_dispersion_tensor(prob.plasma, ω, k; closure = prob.closure)))
        zs, sat[c] = _slice_zeros(alg, f, (pll, pur))
        if !isnothing(refine)
            # Polish on the deflated det (pole-free at ω=0, where the raw det
            # defeats Muller); runaways past the padded window are dropped.
            zs = _polish_seeds(zs, f, refine, 1.0e-4 * diag)
            filter!(ω -> _in_box((pll, pur), ω), zs)
        end
        filter!(ω -> abs(ω) > gate0, zs)
        zv[c], nev[c] = zs, n[]
        return nothing
    end
    nevals = sum(nev)
    stats() = SolveStats(nevals, (time_ns() - t0) / 1.0e9)
    if m == 0
        k0 = kf()
        roots = [
            DispersionBranch(ω, k0, residual(prob, ω, k0))
                for ω in zv[] if _in_box(prob.region, ω)
        ]
        return SurveySolution(roots, stats(), _retcode(roots, any(sat)), prob, alg)
    end
    sheets = link_sheets(grids, zv; gate = @something(cache.gate, diag / 8), maxgap)
    filter!(sh -> any(_in_box(prob.region, ω) for (_, ω) in sh), sheets)
    T = _realtype(prob)
    roots = _tmap(sheets) do sh
        ωs = [Complex{T}(ω) for (_, ω) in sh]
        ks = [kf(p...) for (p, _) in sh]
        res = [residual(prob.plasma, ω, k; closure = prob.closure) for (ω, k) in zip(ωs, ks)]
        DispersionBranch(ωs, ks, res)
    end
    return SurveySolution(roots, stats(), _retcode(roots, any(sat)), prob, alg)
end

_retcode(roots, saturated) =
    isempty(roots) ? :Failure : (saturated ? :Partial : :Success)


# Muller-polish seeds to machine accuracy on `f`, dedup by `dedup`.
# A diverged polish drops its seed — the artifact filter: fit poles with no
# zero of f nearby (Froissart) diverge; every observed genuine stall is the
# structural ω=0 double zero the origin gate removes anyway.
function _polish_seeds(seeds, f, alg::Muller, dedup)
    CT = complex(float(eltype(seeds)))
    out = CT[]
    for ωc in seeds
        h = _seed_offset(ωc)
        ω = muller(f, ωc - h, ωc + h, ωc + im * h; alg.atol, alg.maxiter)
        isfinite(ω) || continue
        any(o -> abs(ω - o) < dedup, out) || push!(out, ω)
    end
    return out
end
