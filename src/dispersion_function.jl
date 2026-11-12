# Fixed-k evaluation layer: per-species susceptibility plans + the callable
# dispersion function root finders consume.
# Everything ω-independent (Γ/Bessel tables, harmonic caps, reduced moments)
# is hoisted to plan time for efficient ω-evaluations.

"""
    DispersionFunction(plasma, k; deflate=true, scale=1.0, mode=:det, backend=FixedNodeEval())

Callable fixed-`k` dispersion function `D(ω) = mode(𝒟)/scale`, zero at the dispersion relation's roots.
Construction plans every species ([`plan_contribution`](@ref)) once for this `k`.

- `deflate=true`: `𝒟` is the deflated `ω̃²𝒟`, pole-free at `ω = 0` (survey form); `false`: the raw `𝒟` (seeded form).
- `mode`: the [`TensorReduction`](@ref).
"""
struct DispersionFunction{L,P,K,M,S,R<:TensorReduction}
    plans::L
    plasma::P
    k::K
    curl::M       # k̃k̃ᵀ − k̃²I
    scale::S
    deflate::Bool
    mode::R
end

DispersionFunction(::Plasma, k; kw...) =
    throw(ArgumentError("DispersionFunction evaluates normalized ω: pass `NormalizedPlasma(plasma)` and `k` in `Ω_ref/c`"))

function DispersionFunction(plasma, k; deflate=true, scale=1.0, mode=FullDet(), backend=FixedNodeEval())
    m = check_reduction(mode, plasma, k)
    plans = map(s -> plan_contribution(s, k; backend), NormalizedPlasma(plasma).populations)
    return DispersionFunction(plans, plasma, k, _curlcurl(k), scale, deflate, m)
end


"""
    plan_contribution(s, k; backend=FixedNodeEval(), kw...) -> plan

Fixed-`k` susceptibility plan: `plan(ω) ≈ contribution(s, ω, k; closure)`, the closure of `backend`.
VDFs override `plan_contribution(vdf, s, k; kw...)` to hoist ω-independent
setup; the fallback closes over the direct path. `scaled_contribution(plan, ω)`
is `ω̃²·plan(ω)`, overridden where the ω→0 poles cancel analytically (cold).
"""
plan_contribution(s, k; kw...) = (s = _solver_population(s); plan_contribution(s.vdf, s, k; kw...))
plan_contribution(vdf, s, k; backend=FixedNodeEval(), kw...) = GenericKPlan(s, k, _closure(backend))

struct GenericKPlan{S,K,C}
    s::S
    k::K
    closure::C
end
(p::GenericKPlan)(ω) = contribution(p.s, ω, p.k; closure=p.closure)
scaled_contribution(p::GenericKPlan, ω) = scaled_contribution(p.s, ω, p.k; closure=p.closure)
scaled_contribution(p, ω) = complex(ω)^2 * p(ω)

function _tensor(D::DispersionFunction, ω)
    ω̃² = complex(ω)^2
    return if D.deflate
        ω̃² * I + _guarded_sum(p -> scaled_contribution(p, ω), D.plans) + D.curl
    else
        I + _guarded_sum(p -> p(ω), D.plans) + D.curl / ω̃²
    end
end
(D::DispersionFunction)(ω) = D.mode(_tensor(D, ω), D.k) / D.scale
residual(D::DispersionFunction, ω) = isfinite(ω) ? _resid(_tensor(D, ω)) : NaN

_hadamard(D) = prod(norm(D[i, :]) for i in 1:3)

# Seeded solvers' raw-tensor function; reduction and scale are fixed at the seed `ω0`.
function seeded_function(plasma, k, ω0, opts)
    D = DispersionFunction(plasma, k; opts..., deflate=false)
    T = _tensor(D, ω0)
    mode = _seeded_mode(D.mode, T)
    s = _hadamard(T)
    scale = isfinite(s) && s > 0 ? s^_scalepow(mode) : one(s)
    return DispersionFunction(D.plans, D.plasma, D.k, D.curl, scale, D.deflate, mode)
end

# `ImplicitDet` only when the seed's least singular value singles out one mode (σ₃/σ₂ < 0.1).
_seeded_mode(m::TensorReduction, T) = m
function _seeded_mode(m::FullDet, T)
    all(isfinite, T) || return m
    # 𝒟 = Σᵢ σᵢuᵢvᵢᴴ ⇒ adj 𝒟 = det 𝒟 · Σᵢ vᵢuᵢᴴ/σᵢ ≈ σ₁σ₂ v₃u₃ᴴ (up to phase), and |det|·‖𝒟‖/‖adj‖² ≈ σ₃/σ₂.
    Ts = T / maximum(abs, T)
    A = _adjugate(Ts)
    abs(det(Ts)) * norm(Ts) < 0.1 * norm(A)^2 || return m
    x = A[:, argmax(j -> norm(A[:, j]), 1:3)]          # ∝ v₃, right null vector
    y = conj(A[argmax(i -> norm(A[i, :]), 1:3), :])    # ∝ u₃, left null vector
    return ImplicitDet(normalize(x), normalize(y))
end

# One root from a guess by `findroot(alg, …; kw...)`, the plasma prepared once for every wavevector.
struct SeedSolver{P, A, O, K}
    plasma::P
    alg::A
    opts::O
    kw::K
end
_seedsolver(prob, alg, opts, kw = (;)) = SeedSolver(_prepare(prob, opts).plasma, alg, opts, kw)

function (s::SeedSolver)(k, ω0)
    D = seeded_function(s.plasma, k, ω0, s.opts)
    ω, nevals, retcode = findroot(s.alg, D, ω0; s.kw...)
    return (; omega = ω, resid = residual(D, ω), nevals, retcode)
end

DispersionFunction(prob::DispersionProblem{<:Region}, k, opts) = DispersionFunction(prob.plasma, k; opts...)

