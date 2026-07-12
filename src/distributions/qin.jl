# --- Evaluator A: Qin closed-orbit  (complex-order Bessel) -----------
# See derivation.md §3A

# Empirically A is a CROSS-VALIDATION backend, not a speedup.
# Using residue extraction so the first integrand is smooth in 2-D (near-resonance peaks removed) and
# the second is a 1-D p⊥ integral carrying the analytic pole + Landau residue
function _coupled_contribution(::Newberger, ::NonRelativistic, c::PreparedVDF, s, ω, k; rtol = 1.0e-6, norm = NORM)
    d = c.vdf
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    lo, hi = d.para
    qlo, qhi = d.perp
    # kz=0: a=ω/Ω is u-independent — no p∥ resonance poles, nothing to peel
    ns = iszero(kz) ? (1:0) : _resonance_harmonics(ω, Ω, kz, lo, hi)
    ζs = [(ω - n * Ω) / kz for n in ns]
    ε = max(qlo, sqrt(eps(real(ω))) * qhi)   # perp lower bound (ε edge-removes the p⊥=0 origin)
    lpole_terms = _lpole_term.(ζs, lo, hi, sign(kz), true)  # u-integral of the analytic pole term
    ρs = Vector{SVector{6, typeof(ω)}}(undef, length(ns))
    β, a0, ak = kperp / Ω, ω / Ω, kz / Ω
    kzω, Ωω = kz / ω, Ω / ω
    residue_prefactor = -2π * Ω / kz
    χ = QuadGK.quadgk(ε, qhi; rtol, norm) do w
        _qin_residues!(ρs, d, ns, ζs, w, β, kzω, residue_prefactor)
        # smooth p∥ remainder: full resummed integrand minus the peeled poles.
        inner = QuadGK.quadgk(lo, hi; rtol, norm) do u
            val = _qin_integrand_nonrel(u, w, d, a0, ak, β, kzω, Ωω)
            @inbounds for i in eachindex(ζs)
                val = val .- ρs[i] .* safe_inv(u - ζs[i])
            end
            val
        end[1]
        @inbounds for i in eachindex(ζs)
            inner = inner .+ ρs[i] .* lpole_terms[i]
        end
        inner
    end[1]
    return (s.Pi2 / (ω * Ω)) .* _antisymmat(χ)
end

# (γ,p∥) edge-mapped relativistic Newberger backend, cross-validation only.
# Valid for Im ω ≥ 0 (and damped ω with no resonance in support); it has no
# damped-side continuation (rim branch-cut / apex).
function _coupled_contribution(::Newberger, ::Relativistic, c::PreparedVDF, s, ω, k; norm = NORM)
    d = c.vdf
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    β = kperp / Ω
    pmax = maximum(abs, d.para)
    γmax = sqrt(1 + pmax^2 + d.perp[2]^2)
    umaxmax = sqrt(γmax^2 - 1)
    # Generous harmonic window for (γ,u)∈[1,γmax]×[−umax,umax] —
    # min/max keeps it valid for either sign of kz, Ω, or Re ω; the per-γ in-range
    # filter discards non-crossing n.
    rs = extrema(
        (real(ω) * γ - kz * u) / Ω for γ in (one(γmax), γmax), u in (-umaxmax, umaxmax)
    )
    nlo, nhi = floor(Int, rs[1]) - 1, ceil(Int, rs[2]) + 1
    σ = sign(kz)
    function inner(γ)
        umax = sqrt(γ^2 - 1)
        ζρ = NamedTuple{(:ζ, :ρ), Tuple{typeof(ω), SVector{6, typeof(ω)}}}[]
        for n in nlo:nhi
            ζ = (ω * γ - n * Ω) / kz
            -umax < real(ζ) < umax || continue
            w = sqrt(complex(γ^2 - 1 - ζ^2))
            ρ = ((2π * _U_cov(d, ζ, w, γ, ω, kz)) * (-Ω / kz)) .* _T_n_bare(n, β * w, ζ, w)
            push!(ζρ, (; ζ, ρ))
        end
        fI(u, w) = (2π * _U_cov(d, u, w, γ, ω, kz)) .* _qin_T_bare((ω * γ - kz * u) / Ω, β * w, u, w)
        # Inner edge map (derivation §5.2.2)
        acc = QuadGK.quadgk(-π / 2, π / 2; rtol = 1.0e-7, norm) do θ
            u, w = umax .* sincos(θ)
            val = fI(u, w)
            for p in ζρ
                val = val .- p.ρ ./ (u - p.ζ)
            end
            w .* val
        end[1]
        for p in ζρ
            acc = acc .+ p.ρ .* _lpole_term(p.ζ, -umax, umax, σ, true)
        end
        return acc
    end
    val = QuadGK.quadgk(zero(real(ω)), one(real(ω)); rtol = 1.0e-6, norm) do q
        γ = 1 + (γmax - 1) * q^2
        (2 * (γmax - 1) * q) .* inner(γ)
    end[1]
    # `fI` carries only the resonant 𝒰·𝓣ₙ; add the same pole-free nonresonant 𝒳_B
    bern = _ee33((s.Pi2 / ω^2) * c.cache.bernstein33)
    return (s.Pi2 / (ω^2 * Ω)) .* _antisymmat(val) .+ bern
end

include("qin_sigmas.jl")

@inline function _qin_integrand_nonrel(u, w, d, a0, ak, β, kzω, Ωω)
    a = a0 - ak * u
    z = β * w
    r = u / w
    dfpe, dfpa = d.dgrad(w, u)
    cross = w * dfpa - u * dfpe
    U = dfpe + kzω * cross
    ee = Ωω * r * cross
    bern = @SVector ComplexF64[0, 0, 0, 0, 0, ee]
    return (2π * U) .* _qin_T_bare(a, z, u, w) .+ (2π * w) .* bern
end

# Integer harmonics n whose resonance Re ζ_n=(Re ω−nΩ)/k∥ lies in [lo,hi]. Count is set by k∥·width/Ω — INDEPENDENT of k⊥.
@inline function _resonance_harmonics(ω, Ω, kz, lo, hi)
    r1, r2 = (real(ω) - kz * hi) / Ω, (real(ω) - kz * lo) / Ω
    return floor(Int, min(r1, r2)):ceil(Int, max(r1, r2))
end

# Residue of the bare resummed integrand at p⊥=w for each in-range resonance n,ζ:
# ρ_n = 2π U(ζ,w) 𝓣_n(z,ζ,w) (−Ω/k∥). U is velocity-form.
@inline function _qin_residues!(ρs, d, ns, ζs, w, β, kzω, prefactor)
    z = β * w
    for (i, (n, ζ)) in enumerate(zip(ns, ζs))
        dfpe, dfpa = d.dgrad(w, ζ)
        U = dfpe + kzω * (w * dfpa - ζ * dfpe)
        ρs[i] = (prefactor * U) .* _T_n_bare(n, z, ζ, w)
    end
    return ρs
end

# 𝓣(a,z)=p⊥²·T(a,z): z=k⊥p⊥/Ω₀, w=p⊥, u=p∥. Assembled from the regularized (σ0,σ1,σD,σJ)
@inline function _qin_T_bare(a, z, u, w)
    σ0, σ1, σD, σJ = qin_sigmas(a, z)
    w2 = w^2
    zw = z * w
    xx = a * σ1 * w2
    xy = im * a * σD * w2
    xz = u * σ1 * zw
    yy = σJ * w2
    yz = im * u * σD * zw
    zz = u^2 * σ0
    return SA[xx, xy, xz, yy, yz, zz]
end

# harmonic tensor 𝓣_n = p⊥²·T_n (R≡p⊥·Rₙ, Jw≡p⊥·Jₙ′, Rₙ≡(n/z)Jₙ=½(J_{n−1}+J_{n+1})
@inline function _T_n_bare(n, z, u, w)
    Jm, Jp = besselj(n - 1, z), besselj(n + 1, z)
    Jn = besselj(n, z)
    R = w * (Jm + Jp) / 2          # p⊥·Rₙ, regular
    Jw = w * (Jm - Jp) / 2         # p⊥·Jₙ′
    uJn = u * Jn
    xx, yy, zz = R * R, Jw * Jw, uJn * uJn
    xy = im * R * Jw
    xz = uJn * R
    zy = im * uJn * Jw
    return SA[xx, xy, xz, yy, zy, zz]
end
