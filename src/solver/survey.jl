# Generic global survey: per-slice zero finding.

include("linking.jl")

struct SurveyCache{P, A, R, L}
    prob::P
    alg::A
    refine::R
    linking::L
end


function _init(
        prob::DispersionProblem{<:Region}, alg;
        refine = Muller(), linking = (;)
    )
    linking = (; gate = _boxdiag(prob.target) / 8, linking...)
    return SurveyCache(prepare(prob), alg, refine, linking)
end

function CommonSolve.solve!(cache::SurveyCache)
    t0 = time_ns()
    (; prob, alg, refine) = cache
    grids = paramgrids(prob.k)
    kf = wavefun(prob.k)
    ks = map(c -> kf(map(getindex, grids, Tuple(c))...), CartesianIndices(map(length, grids)))
    zv = similar(ks, Vector{ComplexF64})
    rv = similar(ks, Vector{Float64})
    nev = similar(ks, Int)
    conv = similar(ks, Bool)
    _unblas() do
        Threads.@threads for i in eachindex(ks)
            zv[i], rv[i], nev[i], conv[i] = _pointroots(prob, alg, refine, ks[i])
        end
    end
    stats = SolveStats(sum(nev), (time_ns() - t0) / 1.0e9)
    return build_solution(cache, ks, zv, rv, stats, all(conv))
end

function _unblas(f)
    n = BLAS.get_num_threads()
    n == 1 && return f()
    BLAS.set_num_threads(1)
    try
        return f()
    finally
        BLAS.set_num_threads(n)
    end
end

# All det(𝒟) zeros at one wavevector: discover → trust gate → polish → filter→ dedupe.
function _pointroots(prob, alg, refine, k)
    region = prob.target
    diag = _boxdiag(region)
    gate0 = _in_box(region) ? _origin_gate(alg, diag) : 0.0
    f0 = DispersionFunction(prob.plasma, k; closure = prob.closure, mode = prob.mode, backend = prob.backend)
    # erase only on the ComplexF64 lattice (one probe eval); exotic eltypes pass through
    iscf = f0((region[1] + region[2]) / 2) isa ComplexF64
    f = iscf ? erase_cf(f0) : f0
    trust(z) = _trusted(prob.plasma, z, k)
    zs, n1, converged = discover(alg, f, region; keep = trust)
    filter!(trust, zs)
    zs, n2 = polish!(f, zs, refine)
    filter!(z -> isfinite(z) && _in_box(region, z) && abs(z) > gate0 && trust(z), zs)
    zs = consolidate(zs; atol = 1.0e-4 * diag)
    return zs, map(z -> residual(f0, z), zs), n1 + n2, converged
end

_trusted(plasma, ω, k) = all(s -> trusted(s.vdf, s, ω, k), NormalizedPlasma(plasma).populations)

# m=0 sweeps flow through as 0-dim arrays; collapse to scalars at construction
_scalarize(x::AbstractArray{<:Any, 0}) = x[]
_scalarize(x) = x

function build_solution(cache::SurveyCache, ks, roots, resids, stats, converged = true)
    (; prob, alg) = cache
    sheets = link(roots; cache.linking...)
    branches = _branches(roots, resids, sheets, _scalarize(ks), _realtype(prob))
    filter!(b -> any(_in_box.(Ref(prob.target), b.omega)), branches)
    return SurveySolution(branches, stats, _retcode(branches, converged), prob, alg)
end

# Inference-friendly by construction: the ::Type{T} barrier keeps T static
function _branches(roots, resids, sheets, ks, ::Type{T}) where {T}
    return map(sheets) do sheet
        ωs = _scalarize(_gather(roots, sheet, Complex{T}(NaN, NaN)))
        DispersionBranch(ωs, ks, _scalarize(_gather(resids, sheet, T(NaN))), true)
    end
end

_gather(vals, sheet, nan) = map((v, i) -> iszero(i) ? nan : oftype(nan, v[i]), vals, sheet)

_boxdiag((ll, ur)) = abs(ur - ll)

polish!(f, ωs, ::Nothing) = ωs, 0

# `converged=false` in the discovery process may miss some roots ⇒ `Saturated`.
_retcode(roots, converged) =
    isempty(roots) ? ReturnCode.Failure :
    converged ? ReturnCode.Success : ReturnCode.Saturated

# Deduplicate candidate roots using absolute distance `atol`
function consolidate(points; atol)
    out = empty(points)
    for point in points
        any(other -> abs(point - other) < atol, out) || push!(out, point)
    end
    return out
end
