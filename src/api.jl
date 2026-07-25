"""
    trusted(vdf, s, ω, k)::Bool

Whether the VDF can resolve a candidate root at `ω`; surveys drop the ones it cannot.
`true` for every exact VDF — an approximate one has spurious zeros of its OWN det,
and no test on the det can tell the two apart.
"""
trusted(vdf, s, ω, k) = true

"""
    prepare(x, closure; kw...) -> x′

One-time setup before repeated evaluations.
Precompute (ω,k)-independent quantities held in [`PreparedVDF`](@ref).
"""
prepare(x, closure; kw...) = x

"""
    parallel_even(vdf)::Bool

Whether `f₀` is guaranteed even in `p∥` (no field-aligned drift or asymmetry).
"""
parallel_even(x) = false

function contribution end
function plan_contribution end

"""
    ChiBackend

Selects how a per-`k` susceptibility plan evaluates `χ(ω)`.

- [`FixedNodeEval`](@ref) (default): hoist every ω-independent quantity onto a
  fixed quadrature grid built once per `k`, so ω-evaluation is cheap arithmetic.
- [`AdaptiveEval`](@ref): re-run the adaptive quadrature at every ω.
"""
abstract type ChiBackend end

"""
    FixedNodeEval(; order=12, nprobe=129)

Fixed-node per-`k` plan backend (default). Both knobs size the plan-time grid:

- `order`: Gauss–Legendre order per panel.
- `nprobe`: samples of the marginal proxy `∫|f₀|` that adaptive panel discovery integrates.
"""
Base.@kwdef struct FixedNodeEval <: ChiBackend
    order::Int = 12
    nprobe::Int = 129
end

"""
Per-ω adaptive-quadrature backend. Slower, but re-refines against each ω's
integrand — the fallback when fixed nodes chosen at plan time under-resolve
joint or per-ω structure a marginal proxy can't see.
"""
struct AdaptiveEval <: ChiBackend end


# Contract: `discover(alg, f, region; keep) -> (zeros, nevals, converged)`
function discover end
