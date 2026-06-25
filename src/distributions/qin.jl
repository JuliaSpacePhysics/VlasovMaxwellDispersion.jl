# --- Evaluator A: Qin closed-orbit  (complex-order Bessel) -----------
# `derivation.md` §3A:. χ is then a single 2-D
# momentum cubature of `2π U 𝓣 + 2πp⊥ Bernstein`, with `U` and
# resonances carried by `1/sin(πa)`, `a=(ωγ−k∥p∥)/Ω₀`.
# For `Im ω>0`, the integer-`a` poles sit off the real plane ⇒ plain nested QuadGK.

# Empirically A is a CROSS-VALIDATION backend, not a speedup.
# Using residue extraction so the first integrand is smooth in 2-D (near-resonance peaks removed) and
# the second is a 1-D p⊥ integral carrying the analytic pole + Landau residue
function _coupled_contribution(::Newberger, ::NonRelativistic, d::CoupledVDF, s, ω, k; rtol = 1.0e-7, norm = x -> maximum(abs, x))
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    lo, hi = d.parlo, d.parhi
    ns = _resonance_harmonics(ω, Ω, kz, lo, hi)
    ζs = [(ω - n * Ω) / kz for n in ns]                     # p⊥-independent (nonrel)
    ε = sqrt(eps(real(ω))) * d.perphi
    # 2-D smooth remainder: full resummed integrand minus the peeled poles.
    function smooth2d(x)
        u, w = x[1], x[2]
        val = _qin_integrand(u, w, one(real(ω)), d, ω, Ω, kz, kperp)
        for (ζ, ρ) in zip(ζs, _qin_residues(d, ns, ζs, w, ω, Ω, kz, kperp))
            val = val .- ρ ./ (u - ζ)
        end
        return val
    end
    bulk = first(hcubature(smooth2d, SVector(lo, ε), SVector(hi, d.perphi); rtol, norm))
    # 1-D p⊥ integral of the analytic pole terms (+ Landau for damped modes).
    function poles1d(w)
        acc = zero(SMatrix{3, 3, ComplexF64})
        for (ζ, ρ) in zip(ζs, _qin_residues(d, ns, ζs, w, ω, Ω, kz, kperp))
            acc = acc .+ ρ .* _landau_logfac(ζ, lo, hi)
        end
        return acc
    end
    poles = isempty(ns) ? zero(bulk) : first(QuadGK.quadgk(poles1d, ε, d.perphi; rtol, norm))
    return SMatrix{3, 3, ComplexF64}((s.Pi2 / (ω * Ω)) .* (bulk .+ poles))
end

function _coupled_contribution(::Newberger, ::Relativistic, d::CoupledVDF, s, ω, k; norm = x -> maximum(abs, x))
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    β = kperp / Ω
    pmax = max(abs(d.parlo), abs(d.parhi))
    γmax = sqrt(1 + pmax^2 + d.perphi^2)
    umaxmax = sqrt(γmax^2 - 1)
    # Generous harmonic window; the per-γ in-range filter discards non-crossing n.
    nlo = floor(Int, (real(ω) - kz * umaxmax) / Ω) - 1
    nhi = ceil(Int, (real(ω) * γmax + kz * umaxmax) / Ω) + 1
    function inner(γ)
        umax = sqrt(γ^2 - 1)
        ζρ = NamedTuple[]
        for n in nlo:nhi
            ζ = (ω * γ - n * Ω) / kz
            -umax < real(ζ) < umax || continue
            w = sqrt(complex(γ^2 - 1 - ζ^2))
            ρ = ((2π * _U_cov(d, ζ, w, γ, ω, kz)) * (-Ω / kz)) .* _T_n_bare(n, β * w, ζ, w)
            push!(ζρ, (; ζ, ρ))
        end
        fI(u, w) = (2π * _U_cov(d, u, w, γ, ω, kz)) .* _qin_T_bare((ω * γ - kz * u) / Ω, β * w, u, w)
        # Inner edge map (derivation §5.2.2)
        acc = first(
            QuadGK.quadgk(-π / 2, π / 2; rtol = 1.0e-7, norm) do θ
                u, w = umax .* sincos(θ)
                val = fI(u, w)
                for p in ζρ
                    val = val .- p.ρ ./ (u - p.ζ)
                end
                w .* val
            end
        )
        for p in ζρ
            acc = acc .+ p.ρ .* _landau_logfac(p.ζ, -umax, umax)
        end
        return acc
    end
    val = first(
        QuadGK.quadgk(zero(real(ω)), one(real(ω)); rtol = 1.0e-6, norm) do q
            γ = 1 + (γmax - 1) * q^2
            (2 * (γmax - 1) * q) .* inner(γ)
        end
    )
    # `fI` carries only the resonant 𝒰·𝓣ₙ; add the same pole-free nonresonant 𝒳_B
    bern = _ee33((s.Pi2 / ω^2) * _bernstein_rel(d, γmax))
    return SMatrix{3, 3, ComplexF64}((s.Pi2 / (ω^2 * Ω)) .* val) .+ bern
end

include("qin_sigmas.jl")

@inline function _qin_integrand(u, w, γ, d, ω, Ω0, kz, kperp)
    a = (ω * γ - kz * u) / Ω0
    β = kperp / Ω0
    z = β * w
    r = u / w
    dfpa, dfpe = d.dpar(w, u), d.dperp(w, u)
    cross = w * dfpa - u * dfpe
    U = dfpe + (kz / (ω * γ)) * cross               # velocity-form numerator (§2)
    ee = (Ω0 / (γ * ω)) * r * cross                  # e∥e∥ Bernstein term (§3)
    bern = @SMatrix ComplexF64[0 0 0; 0 0 0; 0 0 ee]
    return (2π * U) .* _qin_T_bare(a, z, u, w) .+ (2π * w) .* bern
end

# Integer harmonics n whose resonance Re ζ_n=(Re ω−nΩ)/k∥ lies in [lo,hi]. Count is set by k∥·width/Ω — INDEPENDENT of k⊥.
@inline function _resonance_harmonics(ω, Ω, kz, lo, hi)
    r1, r2 = (real(ω) - kz * hi) / Ω, (real(ω) - kz * lo) / Ω
    return floor(Int, min(r1, r2)):ceil(Int, max(r1, r2))
end

# Residue of the bare resummed integrand at p⊥=w for each in-range resonance n,ζ:
# ρ_n = 2π U(ζ,w) 𝓣_n(z,ζ,w) (−Ω/k∥). U is velocity-form.
@inline function _qin_residues(d, ns, ζs, w, ω, Ω, kz, kperp)
    β = kperp / Ω
    z = β * w
    return map(ns, ζs) do n, ζ
        dfpa, dfpe = d.dpar(w, ζ), d.dperp(w, ζ)
        U = dfpe + (kz / ω) * (w * dfpa - ζ * dfpe)
        ((2π * U) * (-Ω / kz)) .* _T_n_bare(n, z, ζ, w)
    end
end

# 𝓣(a,z)=p⊥²·T(a,z): z=k⊥p⊥/Ω₀, w=p⊥, u=p∥. Assembled from the regularized (σ0,σ1,σD,σJ)
@inline function _qin_T_bare(a, z, u, w)
    σ0, σ1, σD, σJ = qin_sigmas(a, z)
    w2 = w^2
    zw = z * w
    return @SMatrix ComplexF64[
        a * σ1 * w2 im * a * σD * w2 u * σ1 * zw
        -im * a * σD * w2 σJ * w2 -im * u * σD * zw
        u * σ1 * zw im * u * σD * zw u^2 * σ0
    ]
end

# harmonic tensor 𝓣_n = p⊥²·T_n (R≡p⊥·Rₙ, Jw≡p⊥·Jₙ′, Rₙ≡(n/z)Jₙ=½(J_{n−1}+J_{n+1})
@inline function _T_n_bare(n, z, u, w)
    Jm, Jp1 = besselj(n - 1, z), besselj(n + 1, z)
    Jn = besselj(n, z)
    R = w * (Jm + Jp1) / 2          # p⊥·Rₙ, regular
    Jw = w * (Jm - Jp1) / 2         # p⊥·Jₙ′
    return @SMatrix ComplexF64[
        R * R im * R * Jw u * R * Jn
        -im * R * Jw Jw * Jw -im * u * Jw * Jn
        u * R * Jn im * u * Jw * Jn u * u * Jn^2
    ]
end
