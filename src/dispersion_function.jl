# Fixed-k evaluation layer: per-species susceptibility plans + the callable
# dispersion function root finders consume.
# Everything ω-independent (Γ/Bessel tables, harmonic caps, reduced moments)
# is hoisted to plan time for efficient ω-evaluations.

"""
    DispersionFunction(plasma, k; closure=HarmonicSum(), deflate=true, scale=1.0, mode=:det)
    DispersionFunction(prob::AbstractDispersionProblem)

Callable fixed-`k` dispersion function `D(ω)` whose zeros are the dispersion relation's roots.
Construction plans every species ([`plan_contribution`](@ref)) once for this `k`.

`deflate=true` evaluates the deflated `det(ω̃²𝒟), pole-free at
`ω=0` (the survey form); `deflate=false` the raw `det𝒟` (the seeded-solver form).
Both are divided by `scale` (the seeded path normalizes by the det magnitude
at the seed).

`mode` selects the [`TensorReduction`](@ref) to evaluate (`:det` default).
"""
struct DispersionFunction{L,P,K,C,M,S,R<:TensorReduction}
    plans::L
    plasma::P
    k::K
    closure::C
    curl::M       # k̃k̃ᵀ − k̃²I
    scale::S
    deflate::Bool
    mode::R
end

function DispersionFunction(plasma, k; closure=HarmonicSum(), deflate=true, scale=1.0, mode=FullDet())
    m = check_reduction(mode, plasma, k)
    plans = map(s -> plan_contribution(s, k; closure), NormalizedPlasma(plasma).species)
    return DispersionFunction(plans, plasma, k, closure, _curlcurl(k), scale, deflate, m)
end


"""
    plan_contribution(s, k; closure=HarmonicSum(), kw...) -> plan

Fixed-`k` susceptibility plan: `plan(ω) == contribution(s, ω, k; closure)`.
VDFs override `plan_contribution(vdf, s, k; kw...)` to hoist ω-independent
setup; the fallback closes over the direct path. `scaled_contribution(plan, ω)`
is `ω̃²·plan(ω)`, overridden where the ω→0 poles cancel analytically (cold).
"""
plan_contribution(s, k; kw...) = plan_contribution(s.vdf, s, k; kw...)
plan_contribution(vdf, s, k; closure=HarmonicSum(), kw...) = GenericKPlan(s, k, closure)

struct GenericKPlan{S,K,C}
    s::S
    k::K
    closure::C
end
(p::GenericKPlan)(ω) = contribution(p.s, ω, p.k; closure=p.closure)
scaled_contribution(p::GenericKPlan, ω) = scaled_contribution(p.s, ω, p.k; closure=p.closure)
scaled_contribution(p, ω) = complex(ω)^2 * p(ω)

function (D::DispersionFunction)(ω)
    ω̃² = complex(ω)^2
    M = if D.deflate
        ω̃² * I + _guarded_sum(p -> scaled_contribution(p, ω), D.plans) + D.curl
    else
        I + _guarded_sum(p -> p(ω), D.plans) + D.curl / ω̃²
    end
    return D.mode(M, D.k) / D.scale
end

"`false` for plans that approximate (fit) the susceptibility; gates polish on `exact(D)`."
isexact(x) = true
isexact(D::DispersionFunction) = all(isexact, D.plans)

function exact(D::DispersionFunction)
    (; plasma, k, closure, scale, deflate, mode) = D
    return deflate ?
           (ω -> mode(wave_dispersion_tensor(plasma, ω, k; closure), k) / scale) :
           (ω -> mode(dispersion_tensor(plasma, ω, k; closure), k) / scale)
end

_hadamard(D) = prod(norm(D[i, :]) for i in 1:3)

function DispersionFunction(prob::DispersionProblem{<:Seed})
    s = _hadamard(dispersion_tensor(prob.plasma, prob.target.omega0, prob.k; closure=prob.closure))
    scale = isfinite(s) && s > 0 ? s^_scalepow(prob.mode) : one(s)
    return DispersionFunction(
        prob.plasma, prob.k;
        closure=prob.closure, deflate=false, scale, mode=prob.mode
    )
end
DispersionFunction(prob::DispersionProblem{<:Region}, k) =
    DispersionFunction(prob.plasma, k; closure=prob.closure, mode=prob.mode)

