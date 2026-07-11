dispersion_function(prob::AbstractDispersionProblem) =
    ω -> det(dispersion_tensor(prob.plasma, ω, prob.k; closure = prob.closure))

# det(ω̃²𝒟)=ω̃⁶·det𝒟 is pole-free at ω=0. The genuine light-term pole there
# would otherwise cancel nearby roots' winding
deflated_dispersion_function(prob::AbstractDispersionProblem) =
    ω -> det(wave_dispersion_tensor(prob.plasma, ω, prob.k; closure = prob.closure))

@inline Base.getproperty(prob::DispersionProblem, s::Symbol) =
    s === :f ? _scaled_dispersion_function(prob) : getfield(prob, s)
@inline Base.getproperty(prob::GlobalDispersionProblem, s::Symbol) =
    s === :f ? deflated_dispersion_function(prob) : getfield(prob, s)
Base.propertynames(prob::AbstractDispersionProblem) =
    (fieldnames(typeof(prob))..., :f)

_hadamard(D) = prod(norm(D[i, :]) for i in 1:3)

# σ(adj D) = {σ₁σ₂, σ₁σ₃, σ₂σ₃} ⇒ σ₃/σ₁ = |det D| / (σ_max(adj D)·σ_max(D)).
# σ_max via the normal matrix is cancellation-safe — only the *smallest* eigenvalue of D'D is not.
_smax(A) = sqrt(max(last(eigvals(Hermitian(A' * A))), zero(real(eltype(A)))))
_adjugate(A::SMatrix{3, 3}) = @inbounds SMatrix{3, 3}(
    A[2, 2] * A[3, 3] - A[2, 3] * A[3, 2], A[2, 3] * A[3, 1] - A[2, 1] * A[3, 3], A[2, 1] * A[3, 2] - A[2, 2] * A[3, 1],
    A[1, 3] * A[3, 2] - A[1, 2] * A[3, 3], A[1, 1] * A[3, 3] - A[1, 3] * A[3, 1], A[1, 2] * A[3, 1] - A[1, 1] * A[3, 2],
    A[1, 2] * A[2, 3] - A[1, 3] * A[2, 2], A[1, 3] * A[2, 1] - A[1, 1] * A[2, 3], A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1],
)
_sigma_ratio(D) = (s = svdvals(D); s[end] / s[1])
function _sigma_ratio(D::SMatrix{3, 3})
    D_scaled = D / maximum(abs, D) # Pre-scaled to avoid overflow in det/adj'adj
    return abs(det(D_scaled)) / (_smax(_adjugate(D_scaled)) * _smax(D_scaled))
end

"""
    residual(plasma, ω, k; closure=HarmonicSum())
    residual(prob, ω)

Relative Eckart–Young distance to singularity, `σ_min(𝒟)/σ_max(𝒟) ∈ [0,1]`,
is used as the scale-invariant residual. `NaN` for non-finite `ω` or `𝒟`.

Known blind spot: as `ω → 0` transverse terms inflate, causing `σ_max ∝ 1/ω²`.
The determinant's structural origin zero reads small but may not be a true root.
"""
function residual(plasma, ω, k; closure = HarmonicSum())
    isfinite(ω) || return NaN
    D = dispersion_tensor(plasma, ω, k; closure)
    all(isfinite, D) || return NaN
    return _sigma_ratio(D)
end
residual(prob::AbstractDispersionProblem, ω, k = prob.k) =
    residual(prob.plasma, ω, k; closure = prob.closure)

function _scaled_dispersion_function(prob::DispersionProblem)
    f = dispersion_function(prob)
    s = _hadamard(dispersion_tensor(prob.plasma, prob.omega0, prob.k; closure = prob.closure))
    return isfinite(s) && s > 0 ? (ω -> f(ω) / s) : f
end

include("solver/muller.jl")
include("solver/GRPF.jl")
include("solver/ArcLength.jl")
include("solver/AAA.jl")
include("solver/survey.jl")

"""
    solve(prob::DispersionProblem, alg = Muller()) -> DispersionSolution
"""
CommonSolve.solve(prob::DispersionProblem{<:Any, <:Wavenumber}) = CommonSolve.solve(prob, Muller())
CommonSolve.solve(prob::DispersionProblem) = CommonSolve.solve(prob, ArcLength())
CommonSolve.init(prob::DispersionProblem, alg; kwargs...) =
    CommonSolve.init(prob, ArcLength(base = alg); kwargs...)


"""
    wave_dispersion_tensor(plasma, ω, k::Wavenumber; closure=HarmonicSum())

Deflated form `ω̃²·𝒟 = ω̃²ε + (k̃k̃ᵀ − k̃²I)`, built as `ω̃²I + ω̃²χ + curlcurl` so 
the light-term `curlcurl/ω̃²` pole for original `det𝒟` and any `χ` pole at `ω=0` 
(cold `ε`'s `1/ω²`, `1/ω` terms) cancel analytically.

Otherwise its winding partially cancels nearby roots, causing GRPF to miss them
and report a spurious net pole.
"""
function wave_dispersion_tensor(plasma, ω, k; kwargs...)
    ω2χ = _guarded_sum(s -> scaled_contribution(s, ω, k; kwargs...), plasma)
    return ω^2 * I + ω2χ + _curlcurl(k)
end

"""
    solve(prob::GlobalDispersionProblem, alg=AAA(); refine=Muller(), kw...)::SurveySolution

Find all root [`DispersionBranch`](@ref)es of the deflated `det(ω̃²𝒟)`: `alg`
([`AAA`](@ref) or [`GRPF`](@ref)) runs at each point of the geometry's
parameter grid, and per-point roots are linked into sheets by
[`link_sheets`](@ref) (`gate` defaults to ⅛ of the box diagonal). 
The ω box is a soft window tracked `pad` past every edge. 
Fixed `k` gives single-point branches.

`refine` (default [`Muller`](@ref); `nothing` keeps raw fit/mesh roots)
polishes each root and filters out candidates with no nearby zero of the det.
"""
CommonSolve.solve(prob::GlobalDispersionProblem; kwargs...) =
    CommonSolve.solve(prob, AAA(); kwargs...)

function _in_box(region, point = 0)
    ll, ur = region
    return real(ll) ≤ real(point) ≤ real(ur) && imag(ll) ≤ imag(point) ≤ imag(ur)
end
