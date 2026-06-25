# bi-Maxwellian: harmonic sum closed by Z/Z' (parallel) and Γ_n=Iₙe^{-λ} (perp).
# Ref: LMV Tensors.jl, Maxwellian specialization

"""
    Maxwellian(; vth_par, vth_perp=vth_par, vd=0)

Drifting bi-Maxwellian, dimensionless thermal speeds `v_th/c` and parallel
drift `vd/c`.
"""
Base.@kwdef struct Maxwellian{T} <: AbstractVDF
    vth_par::T
    vth_perp::T = vth_par
    vd::T = zero(vth_par)
end

Maxwellian(vth_par) = Maxwellian(; vth_par)

@inline thermal_par(d::Maxwellian) = d.vth_par
@inline thermal_perp(d::Maxwellian) = d.vth_perp
@inline drift(d::Maxwellian) = d.vd

function contribution(d::Maxwellian, s, ω, k; rtol=1.0e-8, kwargs...)
    Ω = s.Omega
    kz = para(k)
    kperp = perp(k)
    vthperp = thermal_perp(d)
    ω = complex(float(ω))

    k⊥_Ω = kperp / Ω
    λ = (vthperp^2 / 2) * k⊥_Ω^2

    prefac = s.Pi2 / ω^2
    minharm = 1
    nmax = nmax_bessel(λ)

    f = n -> _maxwellian_harmonic(n, ω, Ω, kz, k⊥_Ω, d.vth_par, vthperp, d.vd)
    χ = converge(f, minharm, rtol; nmax)
    return SMatrix{3,3,ComplexF64}(prefac * χ)
end


# One cyclotron harmonic of the bi-Maxwellian χ. Parallel moments z*F/z*T from Z(ζ); perp moments ⊥* from
# Γ_n=Iₙe^{-λ} and its neighbours via the besselix recurrence.
@inline function _maxwellian_harmonic(n, ω, Ω, kz, k⊥_Ω, vthpar, vthperp, vd)
    nΩ = n * Ω

    # --- parallel integral moments (kz ≠ 0 validated path) ---
    # zᵖF = ∫ vᵖ f∥/(v-ζ); zᵖT = ∫ vᵖ ∂_v f∥/(v-ζ). PDF moments:
    # Z0=Z(ζ), Z1=1+ζZ0, Z2=ζZ1 (∫vⁿe^{-v²}/√π/(v-ζ)).
    σ⁻¹ = 1 / (kz * vthpar)
    ζ = (ω - kz * vd - nΩ) * σ⁻¹
    Z0 = Z(ζ)
    Z1 = 1 + ζ * Z0
    Z2 = ζ * Z1
    a = Z0
    b = Z0 * vd + Z1 * vthpar
    c = Z0 * vd^2 + Z1 * 2 * vthpar * vd + Z2 * vthpar^2
    z0F = -a * σ⁻¹
    z1F = -b * σ⁻¹
    z2F = -c * σ⁻¹
    invth2 = 2 / vthpar^2
    z0T = (z0F * vd - z1F) * invth2
    z1T = (z1F * vd - z2F) * invth2

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
    PF = (-vth²₂) * P∂                 # f⊥′=−v⊥f⊥/vth²₂ ⇒ F slice = −vth²₂·∂F

    z = (z0F, z1F, z2F, z0T, z1T)
    return _chi_mblock(z, P∂, PF, ω, kz, nΩ)
end
