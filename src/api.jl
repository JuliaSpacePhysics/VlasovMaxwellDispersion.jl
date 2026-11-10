"""
    trusted(vdf, s, ω, k)::Bool

Whether the VDF can resolve a candidate root at `ω`; surveys drop the ones it cannot.
`true` for every exact VDF — an approximate one has spurious zeros of its OWN det,
and no test on the det can tell the two apart.
"""
trusted(vdf, s, ω, k) = true

"""
    prepare(x; backend = FixedNodeEval(), kw...) -> x′

One-time setup before repeated evaluations.
Precompute the (ω,k)-independent quantities `backend` ([`ChiBackend`](@ref)) needs, held in [`PreparedVDF`](@ref).
"""
prepare(x; kw...) = x

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
    FixedNodeEval(; order=12, nprobe=129, tol=1e-12)

Fixed-node per-`k` plan backend (default). The knobs size the fixed grid:

- `order`: Gauss–Legendre order per panel.
- `nprobe`: samples of the marginal proxy `∫|f₀|` that adaptive panel discovery integrates.
- `tol`: panels are bisected until the tabulated `∇f₀` interpolates to `tol·max|∇f₀|`; this
  is the accuracy floor of `χ`.

The tables are k-independent and built by `prepare`, which takes the same `backend`.
"""
Base.@kwdef struct FixedNodeEval <: ChiBackend
    order::Int = 12
    nprobe::Int = 129
    tol::Float64 = 1.0e-12
end

"""
    AdaptiveEval(; closure=HarmonicSum())

Per-ω adaptive-quadrature backend. Slower, but re-refines against each ω's
integrand — the fallback when fixed nodes chosen at plan time under-resolve
joint or per-ω structure a marginal proxy can't see.

`closure` evaluates the orbit integrals of `CoupledVDF` and `SplineVDF`.
"""
Base.@kwdef struct AdaptiveEval{C} <: ChiBackend
    closure::C = HarmonicSum()
end

_closure(::ChiBackend) = HarmonicSum()
_closure(b::AdaptiveEval) = b.closure


# Contract: `discover(alg, f, region; keep) -> (zeros, nevals, converged)`
function discover end

# Contract of a seeded algorithm: `findroot(alg, D, ω0; kw...) -> (ω, nevals, retcode)`, `NaN` when the search fails.
function findroot end
