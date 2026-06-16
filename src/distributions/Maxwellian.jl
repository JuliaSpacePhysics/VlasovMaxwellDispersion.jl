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

Continuation(::Maxwellian) = Analytic()

@inline thermal_par(d::Maxwellian) = d.vth_par
@inline thermal_perp(d::Maxwellian) = d.vth_perp
@inline drift(d::Maxwellian) = d.vd

function contribution(d::Maxwellian, s::Species, ω, k; kwargs...)
    Ω = s.Omega
    kz = para(k)
    kperp = perp(k)
    vthpar = thermal_par(d)
    vthperp = thermal_perp(d)
    vd = drift(d)
    ω = complex(float(ω))

    k⊥_Ω = kperp / Ω
    λ = (vthperp^2 / 2) * k⊥_Ω^2

    prefac = s.Pi2 / ω^2
    minharm = 1
    rtol = 1.0e-8
    nmax = nmax_bessel(λ)

    f = n -> _maxwellian_harmonic(n, ω, Ω, kz, kperp, k⊥_Ω, vthpar, vthperp, vd)
    χ = converge(f, minharm, rtol; nmax)
    return SMatrix{3,3,ComplexF64}(prefac * χ)
end


# One cyclotron harmonic of the bi-Maxwellian χ. Parallel moments z*F/z*T from Z(ζ); perp moments ⊥* from
# Γ_n=Iₙe^{-λ} and its neighbours via the besselix recurrence.
@inline function _maxwellian_harmonic(n, ω, Ω, kz, kperp, k⊥_Ω, vthpar, vthperp, vd)
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

    # --- perpendicular integral moments ---
    vth²₂ = vthperp^2 / 2
    λ = vth²₂ * k⊥_Ω^2
    ∫n₋ = Gamma_n(n - 1, λ)
    ∫n₊ = n == 0 ? ∫n₋ : Gamma_n(n + 1, λ)
    ∫n₀ = n == 0 ? Gamma_n(n, λ) : λ / (2n) * (∫n₋ - ∫n₊)

    JF = ∫n₀                                     # ∫Jn² F⊥ 2π v⊥
    J∂F = -∫n₀ / vth²₂                            # ∫Jn² ∂F⊥ 2π
    JdJ∂F = -k⊥_Ω * ((∫n₋ + ∫n₊) / 2 - ∫n₀)       # ∫Jn ∂Jn ∂F⊥ 2π v⊥
    JdJF = -JdJ∂F * vth²₂                         # ∫Jn ∂Jn F⊥ 2π v⊥²
    ∂J²∂F = λ * (∫n₊ - 2∫n₀ + ∫n₋) + n * (∫n₊ - ∫n₋) / 2  # ∫∂Jn² ∂F⊥ 2π v⊥²
    ∂J²F = -vth²₂ * ∂J²∂F                         # ∫∂Jn² F⊥ 2π v⊥³

    z = (z0F, z1F, z2F, z0T, z1T)
    iszero(kperp) &&
        return _chi_mblock_kperp0(z, n, vth²₂, JF, JdJF, JdJ∂F, ∂J²F, ∂J²∂F, ω, kz)
    p = (; JF, J∂F, JdJF, JdJ∂F, ∂J²F, ∂J²∂F)
    return _chi_mblock(z, p, ω, kz, kperp, n / k⊥_Ω)
end


# Parallel propagation k⊥=0 limit of `_chi_mblock` (Maxwellian-only).
@inline function _chi_mblock_kperp0(z, n, vth²₂, JF, JdJF, JdJ∂F, ∂J²F, ∂J²∂F, ω, kz)
    z0F, z1F, z2F, z0T, z1T = z
    D(X, Y, a, b, c) = kz * (-X * a + Y * b) + ω * X * c
    o = abs(n) == 1
    m11 = (o / 2) * (kz * (z1F + vth²₂ * z0T) - ω * z0F)
    m21 = -im * n * m11
    m22 = D(∂J²∂F, ∂J²F, z1F, z0T, z0F)
    m32 = im * D(JdJ∂F, JdJF, z2F, z1T, z1F)
    m33 = ω * JF * z1T
    return @SMatrix ComplexF64[m11 -m21 0; m21 m22 -m32; 0 m32 m33]
end
