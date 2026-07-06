# Kappa-family. Parallel Cauchy moments of (b+u²/a)^{-M} close analytically.
"""
    Kappa(vth, kappa)

1-D kappa factor `∝ (1 + v²/a)^{-(κ+1)}`, `a = (κ-3/2)·vth²`, `κ>3/2`
(→ [`Gaussian`](@ref)`(vth)` as `κ→∞`).
"""
struct Kappa{T, K}
    a::T
    kappa::K

    function Kappa(vth::T, κ::K) where {T, K}
        κ > 1.5 || throw(ArgumentError("Kappa needs κ > 3/2 (finite ⟨v²⟩)"))
        a = (κ - 1.5) * vth^2
        return new{T, K}(a, κ)
    end
end


@inline (p::Kappa)(v) = (1 + v^2 / p.a)^(-(p.kappa + 1))

"""
    ProductBiKappa(; vth_para, vth_perp=vth_para, kappa_para, kappa_perp=kappa_para)

Separable product bi-kappa `f₀ = f⊥(p⊥)·f∥(p∥)` with *independent* spectral indices:

    f₀ ∝ Kappa(vth_perp, kappa_perp) ⊗ Kappa(vth_para, kappa_para)

Note `ProductBiKappa(κ,κ) ≠ BiKappa(κ)`.
"""
ProductBiKappa(; vth_para, vth_perp = vth_para, kappa_para, kappa_perp = kappa_para) =
    Kappa(vth_perp, kappa_perp) ⊗ Kappa(vth_para, kappa_para)


# Base integrals H_m(ζ)=∫p^m/((p²+β²)^M(p-ζ))dp, m=0,1,2. Integer M: order-M pole at p=iβ,
# closed UHP (integrand ~p^{m-2M-1}); the 2πi·ζᵐ/(ζ²+β²)^M pole term is the uniform
# Landau-causal continuation for σ=sign(k∥)>0. For σ<0 the causal side is Im ζ<0, where
# the UHP closure holds only Res_{iβ}.
function _kappa_Hm(ζ, β2, M::Integer, σ = 1)
    iβ = im * sqrt(β2)
    T = promote_type(typeof(ζ), typeof(β2))
    twoiβ, dζ = 2iβ, iβ - ζ
    # Res_{p=iβ} via c_p=[tᵖ](2iβ+t)^{-M}(dζ+t)^{-1}, t=p-iβ. The (dζ+t)^{-1} series is
    # geometric, so the Taylor convolution collapses to prefix sums of A_k=[tᵏ](2iβ+t)^{-M}
    # (ratio recurrence — no binomial overflow):  c_p = (1/dζ)(-1/dζ)ᵖ Σ_{k≤p} A_k(-dζ)ᵏ.
    t = twoiβ^(-M)
    S = t
    S3 = S2 = S                              # S_{M-3}, S_{M-2}; init covers M=3 (p=0)
    for k in 1:(M - 1)
        t *= (M + k - 1) * dζ / (k * twoiβ)
        S += t
        k == M - 3 && (S3 = S)
        k == M - 2 && (S2 = S)
    end
    r = -1 / dζ
    cM1 = S * r^(M - 1) / dζ
    cM2 = S2 * r^(M - 2) / dζ
    cM3 = S3 * r^(M - 3) / dζ
    invden = σ > 0 ? 1 / (ζ^2 + β2)^M : zero(T)   # 2πi·ζᵐ·invden = Landau residue
    pref = 2π * im
    H0 = pref * (cM1 + invden)
    H1 = pref * ((iβ * cM1 + cM2) + ζ * invden)
    H2 = pref * ((2iβ * cM2 + cM3 - β2 * cM1) + ζ^2 * invden)
    return H0, H1, H2
end

# Non-integer M: residue fails at the branch point.
# H₀ is the Mace–Hellberg kappa-Z, a single Gauss ₂F₁ (Euler integral G&R 3.259.3)
function _kappa_Hm(ζ, β2, M, σ = 1)
    H₀ = _kappa_H0(ζ, β2, M, σ)
    N₀ = sqrt(π) * gamma(M - 0.5) / gamma(M) * β2^((1 - 2M) / 2)
    H₁ = N₀ + ζ * H₀
    H₂ = ζ * H₁
    return H₀, H₁, H₂
end

# The library's principal ₂F₁ is the direct integral for Im ζ≥0; the direct value for
# Im ζ<0 is its Schwarz reflection (real kernel). The causal side is σ·Im ζ>0; on the
# Landau-crossed side continue with the jump σ·2πi·g(ζ), g = (ζ²+β²)^{-M}.
function _kappa_H0(ζ, β2, M, σ = 1)
    direct = if imag(ζ) >= 0 && !(σ < 0 && imag(ζ) == 0)   # at real ζ take the σ-home limit
        im * sqrt(oftype(β2, π)) * gamma(M + 0.5) / gamma(M + 1) * β2^(-M) *
            _₂F₁(M, 0.5, M + 1, 1 + ζ^2 / β2)
    else
        conj(_kappa_H0(complex(real(ζ), abs(imag(ζ))), β2, M))
    end
    return σ * imag(ζ) < 0 ? direct + σ * 2π * im / (ζ^2 + β2)^M : direct
end

# M_F=(-1/kz)𝒞[uᵐf∥], M_T=(-1/kz)𝒞[uᵐf∥′]. f∥=C(1+u²/a)^{-(κ+1)} at exponent κ+1;
# f∥′=-2C(κ+1)u/a·(…)^{-(κ+2)} raises it to κ+2
function para_moments(p::Kappa, ω, kz, nΩ)
    κ, a = p.kappa, p.a
    if iszero(kz)
        # no u-pole: (1, 0, ⟨u²⟩=a/(2κ−1), 0, ∫uf′=−1)/Δ
        invΔ = 1 / (complex(ω) - nΩ)
        return (invΔ, zero(invΔ), a / (2κ - 1) * invΔ, zero(invΔ), -invΔ)
    end
    ζ = (ω - nΩ) / kz
    σ = sign(kz)
    H0, H1, H2 = _kappa_Hm(ζ, a, κ + 1, σ)
    _, G1, G2 = _kappa_Hm(ζ, a, κ + 2, σ)
    C = gamma(κ + 1) / (sqrt(π * a) * gamma(κ + 0.5)) * a^(κ + 1)   # 1-D norm × (p²+a)^{-M} rescale
    pf = -C / kz
    tf = 2 * C * (κ + 1) / kz
    return (pf * H0, pf * H1, pf * H2, tf * G1, tf * G2)
end

nmax_harm(p::Kappa, β) = nmax_bessel(β^2 * p.a / (p.kappa - 1) / 2)   # ⟨p⊥²⟩=a/(κ-1)

function perp_moments(p::Kappa, n, β; rtol = 1.0e-8)
    κ, a = p.kappa, p.a
    C = κ / (π * a)                          # 2-D normalization
    vc = sqrt(a / (κ - 1))
    P = QuadGK.quadgk(zero(vc), oftype(vc, Inf); rtol) do v
        K = _perp_Bessel_bilinear(n, β, v)
        D = 1 + v^2 / a
        df = -2C * (κ + 1) * v / a * D^(-(κ + 2))
        vcat(df .* K, (v * C * D^(-(κ + 1))) .* K)
    end[1]
    # bilinear order (11,12,22,13,23,33) → `_symmat` order (11,12,13,22,23,33)
    P∂ = 2π .* _symmat(P[1], P[2], P[4], P[3], P[5], P[6])
    PF = 2π .* _symmat(P[7], P[8], P[10], P[9], P[11], P[12])
    return P∂, PF
end

# Slice moments Gₘ=∫uᵐ(b+u²/a∥)^{-M}/(u-ζ)du = a∥^M·H_m(ζ, a∥b, M), M=κ+2.
@inline function _kappa_Gm(ζ, a_para, b, κ, σ = 1)
    M = κ + 2
    H0, H1, H2 = _kappa_Hm(ζ, a_para * b, M, σ)
    aM = a_para^M
    return aM * H0, aM * H1, aM * H2
end

# kz=0 slice: pole-free plain moments Sₘ=∫uᵐ(b+u²/a∥)^{-M}du (S₁=0 by parity), M=κ+2.
@inline function _kappa_Gm0(a_para, b, κ)
    M = κ + 2
    β2 = a_para * b
    S0 = sqrt(π) * gamma(M - 0.5) / gamma(M) * β2^(0.5 - M)
    S2 = sqrt(π) / 2 * gamma(M - 1.5) / gamma(M) * β2^(1.5 - M)
    aM = a_para^M
    return aM * S0, aM * S2
end
