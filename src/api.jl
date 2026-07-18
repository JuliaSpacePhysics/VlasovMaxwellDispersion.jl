"""
    trusted(vdf, s, ω, k)::Bool

Whether the VDF can resolve a candidate root at `ω`; surveys drop the ones it cannot. `true` for
every exact VDF — an approximate one has exact zeros of its OWN det that are not modes of `f₀`,
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
The perpendicular-propagation factorization (`Ordinary`/`Extraordinary` at `k∥ = 0`) is exact only then: any odd `p∥`
moment couples the `E ∥ B₀` component back to the transverse block.
Defaults to `false` — data-driven VDFs (grids, fits, arbitrary `f₀`) cannot
certify their symmetry; declare a method to opt in.
"""
parallel_even(x) = false

function contribution end
function plan_contribution end


# Contract: `discover(alg, f, region; keep) -> (zeros, nevals, converged)`
function discover end
