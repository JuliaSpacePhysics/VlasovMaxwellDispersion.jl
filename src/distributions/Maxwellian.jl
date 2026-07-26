# Ref: LMV Tensors.jl.

"""
    Maxwellian(; vth, vd=0, vr=nothing)
    Maxwellian(vth)

Drifting bi-Maxwellian (dimensionless speeds `v/c`). 
`vth` is a scalar  (isotropic) or a `(⊥, ∥)` pair. 
A perpendicular ring speed `vr` selects the [`GyroRing`](@ref) `I₀` form.

    f ∝ (Gaussian | GyroRing) ⊗ Gaussian(vth_para, vd)
"""
function Maxwellian(; vth, vd = nothing, vr = nothing)
    vth_perp, vth_para = perp_para(vth)
    perp = isnothing(vr) ? Gaussian(vth_perp) : GyroRing(vth_perp, vr)
    return perp ⊗ Gaussian(vth_para, @something(vd, zero(vth_para)))
end
Maxwellian(vth) = Maxwellian(; vth)

# Perp tensor from the three ring sums Γ_{n-1},Γ_n,Γ_{n+1}. Shared by the direct method
# (Γ via besselix per call — used by the k⊥=0 energy-matched fallback) and the table path.
@inline function _gauss_perp_moments(Γm, Γ0, Γp, vth²₂, λ, n, β)
    Γ′ = (Γm + Γp) / 2 - Γ0
    J∂F = -Γ0 / vth²₂                  # Jn² · ∂F⊥
    JdJ∂F = -β * Γ′                    # Jn Jn′ · ∂F⊥ v⊥
    ∂J²∂F = λ * (Γp - 2Γ0 + Γm) + n * (Γp - Γm) / 2  # Jn′² · ∂F⊥ v⊥²
    RR∂F = n * (Γp - Γm) / 2           # Rn² · ∂F⊥ v⊥²   (= −n²Γ_n/λ)
    RJ∂F = -n * Γ′                     # Rn Jn′ · ∂F⊥ v⊥²
    RnJ∂F = β * (Γp - Γm) / 2          # Rn Jn · ∂F⊥ v⊥
    P∂ = _symmat(RR∂F, RJ∂F, RnJ∂F, ∂J²∂F, JdJ∂F, J∂F)
    PF = (-vth²₂) * P∂
    return P∂, PF
end

# Independent Γ_k(λ) table reuse besselix values across harmonic loop
struct GaussianPerpCtx{T}
    vth²₂::T
    λ::T
    Γ::GammaTable{T}
    nmax::Int
end
function perp_setup(d::Gaussian{<:Any, Nothing}, β)
    vth²₂, β = promote(d.vth^2 / 2, β)
    λ = vth²₂ * β^2
    nmax = nmax_bessel(λ)
    return GaussianPerpCtx(vth²₂, λ, GammaTable(λ, nmax + 1), nmax)  # n+1 reach at outermost harmonic
end
nmax_harm(c::GaussianPerpCtx, β) = c.nmax

@inline function perp_moments(c::GaussianPerpCtx, n, β)
    Γ = c.Γ
    return _gauss_perp_moments(Γ[n - 1], Γ[n], Γ[n + 1], c.vth²₂, c.λ, n, β)
end

"""
    GyroRing(vth, vr)

Gyro-averaged shifted-perp Maxwellian (the `I₀` ring): perp weight is the cold-ring ⊛
Maxwellian convolution `Γ_n^{ring}=Σ_m J_m(Λr)² Γ_{n+m}`, `Λr=k⊥ vr/Ω`.
"""
struct GyroRing{T}
    vth::T
    vr::T
end

GyroRing(vth, vr) = GyroRing(promote(vth, vr)...)

(d::GyroRing)(v) = exp(-((v - d.vr) / d.vth)^2) * besselix(0, 2v * d.vr / d.vth^2)

# Reuse cold-ring spectrum (`J_m(Λr)`-weights) and the Γ_k(λ) table across harmonic loop
struct GyroRingCtx{T}
    σ²::T
    pr::T            # vr
    λ::T
    Λr::T
    Γ::GammaTable{T}
    w0::Vector{T}    # J_m²            (G, G_λ, G_λλ weight)
    w1::Vector{T}    # 2 J_m J_m′      (G_Λ, G_λΛ weight)
    w2::Vector{T}    # 2 J_m′² + 2 J_m J_m″   (G_ΛΛ weight)
    mwin::Int
    nmax::Int
end

# k⊥=0: convolution's n/β factors are singular, but only ⟨v⊥²⟩=vth²+vr² enters χ there →
# energy-matched Gaussian (plain `Gaussian(vth)` would wrongly drop vr²).
function perp_setup(d::GyroRing, β)
    iszero(β) && return perp_setup(Gaussian(sqrt(d.vth^2 + d.vr^2)), β)
    σ², pr, β = promote(d.vth^2 / 2, d.vr, β)
    λ = σ² * β^2
    Λr = β * pr
    mwin = nmax_bessel(Λr^2 / 2)                       # cold-ring convolution window
    nmax = nmax_bessel(λ) + mwin
    Γ = GammaTable(λ, nmax + mwin + 2)                 # k+2 reach at the outermost harmonic
    J(j) = besselj(j, Λr)
    w0 = zeros(typeof(Λr), 2mwin + 1)
    w1, w2 = similar(w0), similar(w0)
    for (i, m) in enumerate(-mwin:mwin)
        Jm = J(m)
        Jmd = (J(m - 1) - J(m + 1)) / 2
        Jmdd = (J(m - 2) - 2J(m) + J(m + 2)) / 4
        w0[i] = Jm^2
        w1[i] = 2 * Jm * Jmd
        w2[i] = 2 * Jmd^2 + 2 * Jm * Jmdd
    end
    return GyroRingCtx(σ², pr, λ, Λr, Γ, w0, w1, w2, mwin, nmax)
end
nmax_harm(c::GyroRingCtx, β) = c.nmax

# Perp moments via the cold-ring ⊛ Maxwellian closure (docs/Maxwellian.md "Ring
# generalization"). Every perp entry reduces to two scalar fundamentals per slice — the
# base Γ_n^{ring} and the v⊥²-moment K — built from one convolution `G` and its (λ,Λr)
# partials, all sharing the precomputed cold-ring spectrum. K = (2σ²+pr²)G + 2σ²λ G_λ + 2σ²Λr G_Λ.
function perp_moments(c::GyroRingCtx, n, β)
    σ², pr, λ, Λr, mwin = c.σ², c.pr, c.λ, c.Λr, c.mwin
    Γ = c.Γ
    Gm = Gz = Gp = Gλ = Gλλ = GΛ = GλΛ = GΛΛ = zero(λ)
    @inbounds for (i, m) in enumerate(-mwin:mwin)
        k = n + m
        Γ0 = Γ[k]
        Γp = Γ[k + 1]
        Γm = Γ[k - 1]
        Γ′ = (Γp + Γm) / 2 - Γ0
        Γp′ = (Γ[k + 2] + Γ0) / 2 - Γp
        Γm′ = (Γ0 + Γ[k - 2]) / 2 - Γm
        Γ′′ = (Γp′ + Γm′) / 2 - Γ′ # Γ_n''(λ)

        w0 = c.w0[i]
        Gm += w0 * Γm
        Gz += w0 * Γ0
        Gp += w0 * Γp
        Gλ += w0 * Γ′
        Gλλ += w0 * Γ′′
        GΛ += c.w1[i] * Γ0
        GλΛ += c.w1[i] * Γ′
        GΛΛ += c.w2[i] * Γ0
    end

    β2 = β^2
    ∂βG = 2σ² * β * Gλ + pr * GΛ
    ∂ββG = 2σ² * Gλ + 4 * σ²^2 * β2 * Gλλ + 4σ² * β * pr * GλΛ + pr^2 * GΛΛ
    K = (2σ² + pr^2) * Gz + 2σ² * λ * Gλ + 2σ² * Λr * GΛ
    Kpar = -K / σ² + pr * β * GΛ + pr^2 * Gz / σ²

    AF = Gz
    BF = ∂βG / 2
    CF = ∂ββG / 2 + ∂βG / (2β) - (n^2 / β2) * Gz + K
    Apar = iszero(n) ? zero(Gz) : (β2 / (2n)) * (Gp - Gm)   # n=0 entry killed by nΩ in zz
    Bpar = -β * Gλ
    ∂βBpar = -Gλ - 2σ² * β2 * Gλλ - β * pr * GλΛ
    Cpar = ∂βBpar + Bpar / β - (n^2 / β2) * Apar + Kpar

    nβ = n / β
    P∂ = _symmat(nβ^2 * Apar, nβ * Bpar, nβ * Apar, Cpar, Bpar, Apar)
    PF = _symmat(nβ^2 * AF, nβ * BF, nβ * AF, CF, BF, AF)
    return P∂, PF
end
