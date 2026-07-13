# Generic global survey: per-slice zero finding + cross-slice linking. Any alg
# implementing `discover(alg, f, region) -> (zeros, saturated)` and
# `_origin_gate(alg, diag)` solves GlobalDispersionProblem at any swept
# dimension through this one path.

include("linking.jl")

struct SurveyCache{P, A, R, L}
    prob::P
    alg::A
    refine::R
    linking::L
end


function CommonSolve.init(
        prob::GlobalDispersionProblem, alg;
        refine = Muller(), linking = (;)
    )
    linking = (; gate = _boxdiag(prob.region) / 8, linking...)
    return SurveyCache(prepare(prob), alg, refine, linking)
end

function CommonSolve.solve!(cache::SurveyCache)
    t0 = time_ns()
    (; prob, alg, refine) = cache
    grids = paramgrids(prob.geometry)
    kf = wavefun(prob.geometry)
    ks = map(c -> kf(map(getindex, grids, Tuple(c))...), CartesianIndices(map(length, grids)))
    zv = similar(ks, Vector{ComplexF64})
    nev = similar(ks, Int)
    Threads.@threads for i in eachindex(ks)
        zv[i], nev[i] = _pointroots(prob, alg, refine, ks[i])
    end
    stats = SolveStats(sum(nev), (time_ns() - t0) / 1.0e9)
    return build_solution(cache, ks, zv, stats)
end


# All det(𝒟) zeros at one wavevector: discover → polish → box/origin gate → dedupe.
function _pointroots(prob, alg, refine, k)
    region = prob.region
    diag = _boxdiag(region)
    gate0 = _in_box(region) ? _origin_gate(alg, diag) : 0.0
    f0 = ω -> det(wave_dispersion_tensor(prob.plasma, ω, k; closure = prob.closure))
    # erase only on the ComplexF64 lattice (one probe eval); exotic eltypes pass through
    f = f0((region[1] + region[2]) / 2) isa ComplexF64 ? erase_cf(f0) : f0
    zs, n1 = discover(alg, f, region)
    zs, n2 = polish!(f, zs, refine)
    filter!(z -> isfinite(z) && _in_box(region, z) && abs(z) > gate0, zs)
    return consolidate(zs; atol = 1.0e-4 * diag), n1 + n2
end

# m=0 sweeps flow through as 0-dim arrays; collapse to scalars at construction
_scalarize(x::AbstractArray{<:Any, 0}) = x[]
_scalarize(x) = x

function build_solution(cache::SurveyCache, ks, values, stats)
    (; prob, alg) = cache
    sheets = link(values; cache.linking...)
    filter!(sh -> any(_in_box.(Ref(prob.region), sh)), sheets)
    roots = _branches(prob, sheets, _scalarize(ks), _realtype(prob))
    return SurveySolution(roots, stats, _retcode(roots), prob, alg)
end

# Inference-friendly by construction: the ::Type{T} barrier keeps T static
function _branches(prob, sheets, ks, ::Type{T}) where {T}
    return map(sheets) do sheet
        ωs = _scalarize(Complex{T}.(sheet))
        DispersionBranch(ωs, ks, _residuals(prob, ωs, ks, T))
    end
end

_residuals(prob, ω::Complex, k, ::Type{T}) where {T} =
    isfinite(ω) ? residual(prob.plasma, ω, k; closure = prob.closure) : T(NaN)
_residuals(prob, ωs, ks, ::Type{T}) where {T} =
    broadcast!((ω, k) -> _residuals(prob, ω, k, T), similar(ωs, T), ωs, ks)

_boxdiag((ll, ur)) = abs(ur - ll)

polish!(f, ωs, ::Nothing) = ωs, 0

_retcode(roots) = isempty(roots) ? :Failure : :Success

# Deduplicate candidate roots using absolute distance `atol`
function consolidate(points; atol)
    out = empty(points)
    for point in points
        any(other -> abs(point - other) < atol, out) || push!(out, point)
    end
    return out
end
