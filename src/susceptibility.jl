"""
    contribution(species/vdf, ω, k)

Susceptibility χ_s(ω,k) from one species or vdf.
"""
@inline contribution(s, ω, k; kwargs...) = contribution(s.vdf, s, ω, k; kwargs...)

function contribution(vdf::AbstractVDF, ω, k; kw...)
    return contribution(NormalizedSpecies(1.0, 1.0, vdf), ω, k; kw...)
end

# ω̃²·χ_s
@inline scaled_contribution(s, ω, k; kwargs...) =
    scaled_contribution(s.vdf, s, ω, k; kwargs...)
scaled_contribution(vdf, s, ω, k; kwargs...) =
    complex(ω)^2 * contribution(vdf, s, ω, k; kwargs...)


# Quadrature-based χ paths raise QuadGK's DomainError when overflow.
# Return NaN so root-finders and ω-scans reject point instead of crashing.
@inline function _guarded_sum(f, plasma)
    _nan_tensor() = SMatrix{3, 3, ComplexF64}(ntuple(_ -> complex(NaN, NaN), 9))
    return try
        mapreduce(f, +, NormalizedPlasma(plasma))
    catch err
        err isa DomainError && isdefined(err, :msg) &&
            startswith(err.msg, "integrand produced") || rethrow()
        _nan_tensor()
    end
end

# ε = I + Σ_s χ_s(ω,k)
dielectric(plasma, ω, k; kwargs...) =
    _guarded_sum(s -> contribution(s, ω, k; kwargs...), plasma) + I

# Curl-curl operator k̃k̃ᵀ - k̃²I
function _curlcurl(k)
    kv = vec3(k)
    return kv * kv' - abs2(k) * I
end

"""
    dispersion_tensor(plasma, ω, k; kw...)

`𝒟(ω,k) = ε + (k̃k̃ᵀ - k̃²I)/ω̃²`. `det(𝒟)=0` is the dispersion relation.
"""
function dispersion_tensor(plasma, ω, k; kwargs...)
    return dielectric(plasma, ω, k; kwargs...) + _curlcurl(k) / complex(ω)^2
end

"""
    wave_dispersion_tensor(plasma, ω, k; kw...)

Deflated form `ω̃²·𝒟 = ω̃²ε + (k̃k̃ᵀ − k̃²I)` so
the light-term `curlcurl/ω̃²` pole and any `χ` pole at `ω=0` 
(cold `ε`'s `1/ω²`, `1/ω` terms) cancel analytically.

Otherwise its winding partially cancels nearby roots =>
coarse-grid root-finding (e.g., GRPF) may miss them and
report a spurious net pole.
"""
function wave_dispersion_tensor(plasma, ω, k; kwargs...)
    ω2χ = _guarded_sum(s -> scaled_contribution(s, ω, k; kwargs...), plasma)
    return ω^2 * I + ω2χ + _curlcurl(k)
end

"Aliases for `dispersion_tensor`"
const 𝒟 = dispersion_tensor

"""
    electrostatic_det(plasma, ω, k::Wavenumber)

Cheap longitudinal path `k̃ · ε · k̃`; its zeros are the electrostatic modes.
"""
function electrostatic_det(plasma, ω, k::Wavenumber; kwargs...)
    ε = dielectric(plasma, ω, k; kwargs...)
    kv = vec3(k)
    return dot(kv, ε, kv)
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
