# Global survey: per-wavevector zero finding; PlasmaDispersionBase links the roots into branches.

const _Discover = Union{AAA, GRPF}

struct PointRoots{P, A, O}
    prob::P
    alg::A
    opts::O
end
function pointsolver(prob, alg::_Discover; kw...)
    isnothing(prob.target.box) &&
        throw(ArgumentError("$(nameof(typeof(alg))) searches a bounded ω box: it solves `(lowerleft, upperright)` targets only"))
    opts = _options(kw)
    return PointRoots(_prepare(prob, opts), alg, opts)
end

# All det(𝒟) zeros at one wavevector: discover → trust gate → polish → filter → dedupe.
function (s::PointRoots)(k)
    (; prob, alg, opts) = s
    region = prob.target
    diag = abs(region[2] - region[1])
    gate0 = 0 in window(region) ? _origin_gate(alg, diag) : 0.0
    f0 = DispersionFunction(prob, k, opts)
    # erase only on the ComplexF64 lattice (one probe eval); exotic eltypes pass through
    iscf = f0((region[1] + region[2]) / 2) isa ComplexF64
    f = iscf ? erase_cf(f0) : f0
    trust(z) = _trusted(prob.plasma, z, k)
    zs, n1, converged = discover(alg, f, region; keep = trust)
    keep(z) = isfinite(z) && abs(z) > gate0 && trust(z)
    filter!(keep, zs)
    zs, n2 = polish!(f, zs, alg.refine)
    filter!(keep, zs)
    zs = consolidate(zs; atol = 1.0e-4 * diag)
    return (; omega = zs, resid = map(z -> residual(f0, z), zs), nevals = n1 + n2, converged)
end

_trusted(plasma, ω, k) = all(s -> trusted(s.vdf, s, ω, k), plasma.populations)

polish!(f, ωs, ::Nothing) = ωs, 0

# Deduplicate candidate roots using absolute distance `atol`
function consolidate(points; atol)
    out = empty(points)
    for point in points
        any(other -> abs(point - other) < atol, out) || push!(out, point)
    end
    return out
end
