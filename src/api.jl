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

function contribution end
function plan_contribution end


# Contract: `discover(alg, f, region; keep) -> (zeros, nevals, converged)`
function discover end
