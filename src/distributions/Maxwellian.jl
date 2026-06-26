# bi-Maxwellian: harmonic sum closed by Z/Z' (parallel) and Γ_n=Iₙe^{-λ} (perp).
# Ref: LMV Tensors.jl, Maxwellian specialization

"""
    Maxwellian(; vth_par, vth_perp=vth_par, vd=0, vr=0)

Drifting bi-Maxwellian, dimensionless thermal speeds `v_th/c` and parallel
drift `vd/c`. A nonzero **perpendicular ring speed** `vr` turns it into a
gyrotropic *ring* (gyro-averaged shifted-perp Maxwellian, `docs/Maxwellian.md`):
the perp weight becomes the cold-ring ⊛ Maxwellian convolution
`Γ_n^{ring}=Σ_m J_m(Λr)² Γ_{n+m}`, `Λr=k⊥ vr/Ω`.
"""
Base.@kwdef struct Maxwellian{T} <: AbstractVDF
    vth_par::T
    vth_perp::T = vth_par
    vd::T = zero(vth_par)
    vr::T = zero(vth_par)
end

Maxwellian(vth_par) = Maxwellian(; vth_par)

@inline thermal_par(d::Maxwellian) = d.vth_par
@inline thermal_perp(d::Maxwellian) = d.vth_perp
@inline drift(d::Maxwellian) = d.vd
@inline ring(d::Maxwellian) = d.vr

function contribution(d::Maxwellian, s, ω, k; rtol = 1.0e-8, kwargs...)
    Ω = s.Omega
    kz = para(k)
    kperp = perp(k)
    vthperp = thermal_perp(d)
    ω = complex(float(ω))

    k⊥_Ω = kperp / Ω
    λ = (vthperp^2 / 2) * k⊥_Ω^2

    prefac = s.Pi2 / ω^2
    minharm = 1

    if iszero(d.vr) || iszero(kperp)
        nmax = nmax_bessel(λ)
        f = n -> _maxwellian_harmonic(n, ω, Ω, kz, k⊥_Ω, d.vth_par, vthperp, d.vd)
    else
        Λr = k⊥_Ω * d.vr
        mwin = nmax_bessel(Λr^2 / 2)                 # cold-ring convolution window
        nmax = nmax_bessel(λ) + mwin
        f = n -> _maxwellian_ring_harmonic(n, ω, Ω, kz, k⊥_Ω, d.vth_par, vthperp, d.vd, d.vr, λ, Λr, mwin)
    end
    χ = converge(f, minharm, rtol; nmax)
    return prefac * χ
end


# One cyclotron harmonic of the bi-Maxwellian χ. Parallel moments M_F/M_T from Z(ζ); perp moments ⊥* from
# Γ_n=Iₙe^{-λ} and its neighbours via the besselix recurrence.
@inline function _maxwellian_harmonic(n, ω, Ω, kz, k⊥_Ω, vthpar, vthperp, vd)
    nΩ = n * Ω

    M = _gaussian_par_moments(ω, kz, nΩ, vthpar, vd)

    # --- perpendicular integral moments (Gaussian closed forms) ---
    # Γ_n ring sums; ring kernel moments (RR/RJ/RnJ, derivation §5.1) use the
    # division-free recurrences (n/λ)Γ_n=(Γ_{n−1}−Γ_{n+1})/2 and Γ_n′=(Γ_{n−1}+Γ_{n+1})/2−Γ_n
    vth²₂ = vthperp^2 / 2
    λ = vth²₂ * k⊥_Ω^2
    Γm = Gamma_n(n - 1, λ)
    Γp = n == 0 ? Γm : Gamma_n(n + 1, λ)
    Γ0 = n == 0 ? Gamma_n(n, λ) : λ / (2n) * (Γm - Γp)
    Γ′ = (Γm + Γp) / 2 - Γ0

    # ∂F slice, indexed as the symmetric outer product of (Rn, Jn′, Jn):
    #   [Rn²  RnJn′  RnJn; ·  Jn′²  Jn′Jn; ·  ·  Jn²]
    J∂F = -Γ0 / vth²₂                  # Jn² · ∂F⊥
    JdJ∂F = -k⊥_Ω * Γ′                 # Jn Jn′ · ∂F⊥ v⊥
    ∂J²∂F = λ * (Γp - 2Γ0 + Γm) + n * (Γp - Γm) / 2  # Jn′² · ∂F⊥ v⊥²
    RR∂F = n * (Γp - Γm) / 2           # Rn² · ∂F⊥ v⊥²   (= −n²Γ_n/λ)
    RJ∂F = -n * Γ′                     # Rn Jn′ · ∂F⊥ v⊥²
    RnJ∂F = k⊥_Ω * (Γp - Γm) / 2       # Rn Jn · ∂F⊥ v⊥
    P∂ = _symmat(RR∂F, RJ∂F, RnJ∂F, ∂J²∂F, JdJ∂F, J∂F)
    PF = (-vth²₂) * P∂

    return _chi_mblock(M, P∂, PF, ω, kz, nΩ)
end


# One harmonic of the gyrotropic *ring* χ (perp ring speed `pr=vr`). Parallel
# moments are identical to the bi-Maxwellian; the perp Bessel moments are the
# cold-ring ⊛ Maxwellian closure (docs/Maxwellian.md "Ring generalization").
#
# Every perp tensor entry reduces to two scalar fundamentals per slice — the base
# Γ_n^{ring} and the v⊥²-moment K — built from one convolution `G` and its (λ,Λr)
# partials, all sharing the cold-ring spectrum `J_m(Λr)²`. The v⊥²-moment closes via
#   K = (2σ²+pr²)G + 2σ²λ G_λ + 2σ²Λr G_Λ ,   σ²=v_thperp²/2.
@inline function _maxwellian_ring_harmonic(n, ω, Ω, kz, β, vthpar, vthperp, vd, pr, λ, Λr, mwin)
    nΩ = n * Ω

    M = _gaussian_par_moments(ω, kz, nΩ, vthpar, vd)

    # --- cold-ring ⊛ Maxwellian perp convolutions (single m-sum) ---
    # G≡Γ_n^{ring}=Σ_m J_m(Λr)² Γ_{n+m}; Gλ/Gλλ replace Γ_{n+m} by Γ'/Γ''; the GΛ
    # family weights by 2J_m J_m′ (the ∂_Λr derivative of J_m²).
    σ² = vthperp^2 / 2
    Gm = Gz = Gp = Gλ = Gλλ = GΛ = GλΛ = GΛΛ = zero(λ)
    for m in -mwin:mwin
        Jm = besselj(m, Λr)
        Jm2 = Jm^2
        Jmd = (besselj(m - 1, Λr) - besselj(m + 1, Λr)) / 2
        Jmdd = (besselj(m - 2, Λr) - 2 * besselj(m, Λr) + besselj(m + 2, Λr)) / 4
        k = n + m
        Γ0 = Gamma_n(k, λ)
        Γp = Gamma_n(k + 1, λ)
        Γpp = Gamma_n(k + 2, λ)
        Γm = Gamma_n(k - 1, λ)
        Γmm = Gamma_n(k - 2, λ)
        Γ′ = (Γp + Γm) / 2 - Γ0
        Γp′ = (Γpp + Γ0) / 2 - Γp
        Γm′ = (Γ0 + Γmm) / 2 - Γm
        Γ′′ = (Γp′ + Γm′) / 2 - Γ′ # Γ_n''(λ)

        Gm += Jm2 * Γm
        Gz += Jm2 * Γ0
        Gp += Jm2 * Γp
        Gλ += Jm2 * Γ′
        Gλλ += Jm2 * Γ′′
        wΛ = 2 * Jm * Jmd
        GΛ += wΛ * Γ0
        GλΛ += wΛ * Γ′
        GΛΛ += (2 * Jmd^2 + 2 * Jm * Jmdd) * Γ0
    end

    β2 = β^2
    ∂βG = 2σ² * β * Gλ + pr * GΛ
    ∂ββG = 2σ² * Gλ + 4 * σ²^2 * β2 * Gλλ + 4σ² * β * pr * GλΛ + pr^2 * GΛΛ
    K = (2σ² + pr^2) * Gz + 2σ² * λ * Gλ + 2σ² * Λr * GΛ
    Kpar = -K / σ² + pr * β * GΛ + pr^2 * Gz / σ²

    # F-slice (×v⊥f⊥) and ∂F-slice (×f⊥′) fundamentals A,B,C (Jₙ², p⊥JₙJₙ′, p⊥²Jₙ′²)
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
    return _chi_mblock(M, P∂, PF, ω, kz, nΩ)
end
