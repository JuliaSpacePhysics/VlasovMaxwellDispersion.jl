# Kappa-family. Parallel Cauchy moments of (b+u²/a)^{-M} close analytically.
"""
    Kappa(θ, kappa) ∝ (1 + v²/(κθ²))^{-(κ+1)}

1-D kappa factor with `a = κ·θ²` → `Gaussian(θ)` as `κ→∞`.
"""
struct Kappa{T, K}
    a::T
    kappa::K

    function Kappa(θ::T, κ::K; check = true) where {T, K}
        check && (κ > 0.5 || throw(ArgumentError("Kappa needs κ > 1/2 (finite ⟨v²⟩)")))
        a = κ * θ^2
        return new{T, K}(a, κ)
    end
end


@inline (p::Kappa)(v) = (1 + v^2 / p.a)^(-(p.kappa + 1))

"""
    ProductBiKappa(; vth_para, vth_perp=vth_para, kappa_para, kappa_perp=kappa_para)

Separable product bi-kappa `f₀ = f⊥(p⊥)·f∥(p∥)` with *independent* spectral indices
and temperature-preserving θ's (`⟨p∥²⟩ = vth_para²/2`, `⟨p⊥²⟩ = vth_perp²`):

    θ∥² = (1 - 1/2κ∥)·vth_para²  (κ∥ > 1/2),   θ⊥² = (1 - 1/κ⊥)·vth_perp²  (κ⊥ > 1)

Note `ProductBiKappa(κ,κ) ≠ BiKappa(κ)`.
"""
function ProductBiKappa(; vth_para, vth_perp = vth_para, kappa_para, kappa_perp = kappa_para)
    kappa_para > 0.5 || throw(ArgumentError("ProductBiKappa needs κ∥ > 1/2 (finite ⟨p∥²⟩)"))
    kappa_perp > 1 || throw(ArgumentError("ProductBiKappa needs κ⊥ > 1 (finite ⟨p⊥²⟩)"))
    return Kappa(sqrt(1 - 1 / kappa_perp) * vth_perp, kappa_perp) ⊗
        Kappa(sqrt(1 - 1 / (2kappa_para)) * vth_para, kappa_para)
end


# Residue-sum assembly shared by every meromorphic parallel factor with a single closed-half-
# plane pole p₀ of order M
# Hₘ = ∫pᵐf/(p−ζ)dp = 2πi·(Res_{p₀}[pᵐf/(p−ζ)] + [σ>0] ζᵐf(ζ)),   m=0,1,2
# σ<0 drops the Landau (ζ-side) residue
@inline function _residue_Hm(p₀, ζ, cM1, cM2, cM3, fζ, σ)
    L = σ > 0 ? fζ : zero(fζ)
    pref = 2π * im
    H0 = pref * (cM1 + L)
    H1 = pref * (p₀ * cM1 + cM2 + ζ * L)
    H2 = pref * (p₀^2 * cM1 + 2p₀ * cM2 + cM3 + ζ^2 * L)
    return H0, H1, H2
end

# Scaled base integrals  aᴹ·Hₘ(ζ) = ∫pᵐ(b+p²/a)^{-M}/(p−ζ)dp,  m=0,1,2  (β²=a·b, pole p=iβ),
# closed UHP. The aᴹ scaling avoids overflow/underflow (aᴹ and β²^{-M} terms) for large κ≳150.
# Valid for integer M ≥ 1.
#
# Residue at iβ: aᴹ·cMⱼ = (−1)^{…}·dζ^{j−1}·Σₖ Dₖ with Dₖ = ρᴹpₖwᵏ, ρ=a/(2iβ·dζ), w=dζ/(2iβ),
# pₖ=Γ(M+k)/(Γ(M)k!). The prefix sums Σ_{k≤M−1}, ≤M−2, ≤M−3 give cM1,cM2,cM3 (Sₚ=0 for p<0 ⇒
# M<3 guards).
function _kappa_Hm_scaled(ζ, a, b, M::Integer, σ = 1)
    β2 = a * b
    iβ = im * sqrt(β2)
    twoiβ, dζ = 2iβ, iβ - ζ
    ρ, w = a / (twoiβ * dζ), dζ / twoiβ
    D = ρ^M

    # Plain forward sweep using ratio recurrence Dₖ=D_{k−1}(M+k−1)w/k from D₀=ρᴹ
    if isfinite(D) && abs(D) >= floatmin(Float64)   # sub-normal/0/Inf D₀ ⇒ anchor instead
        S1 = D
        S2 = M >= 2 ? D : zero(D)           # S_{M-2}=Σ_{k≤M-2}: init S₀
        S3 = M >= 3 ? D : zero(D)           # S_{M-3}
        for k in 1:(M - 1)
            D *= (M + k - 1) * w / k
            S1 += D
            k == M - 2 && (S2 = S1)
            k == M - 3 && (S3 = S1)
        end
        if isfinite(S1)                     # a mid-sweep overflow ⇒ anchor
            s = iseven(M) ? -one(ρ) : one(ρ) # (−1)^{M−1}
            L = (a / (ζ^2 + β2))^M
            return _residue_Hm(iβ, ζ, s * S1, -s * dζ * S2, s * dζ^2 * S3, L, σ)
        end
    end
    # Defer to `_kappa_Hm_anchored`, valid when peak term |Dₖ*| (k*≈|w|(M−1)/(1−|w|)) is still representable
    return _kappa_Hm_anchored(ζ, a, b, M, σ)
end

# Start sweep at k* via loggamma so every partial stays ~O(answer).
function _kappa_Hm_anchored(ζ, a, b, M::Integer, σ = 1)
    β2 = a * b
    iβ = im * sqrt(complex(β2))
    twoiβ, dζ = 2iβ, iβ - ζ
    ρ, w = a / (twoiβ * dζ), dζ / twoiβ
    logρ, logw, lgM = log(ρ), log(w), loggamma(float(M))
    c = M * logρ - lgM                        # log Dₖ = c + loggamma(M+k) − loggamma(k+1) + k·logw
    aw = abs(w)
    kstar = aw >= 1 ? M - 1 : clamp(round(Int, aw * (M - 1) / (1 - aw)), 0, M - 1)
    Dstar = exp(c + loggamma(M + kstar) - loggamma(kstar + 1.0) + kstar * logw)
    tol = 1.0e-300
    SD1 = Dstar
    D = Dstar
    for k in (kstar + 1):(M - 1)
        D *= (M + k - 1) * w / k
        SD1 += D
        abs(D) < tol * abs(SD1) && break
    end
    D = Dstar
    for k in kstar:-1:1
        D *= k / ((M + k - 1) * w)
        SD1 += D
        abs(D) < tol * abs(SD1) && break
    end
    Dm1 = exp(c + loggamma(2M - 1) - lgM + (M - 1) * logw)   # top term k=M-1 (S_{M-2}=SD1−D_{M-1})
    SD2 = SD1 - Dm1
    SD3 = M >= 2 ? SD1 - Dm1 - exp(c + loggamma(2M - 2) - loggamma(M - 1.0) + (M - 2) * logw) : zero(SD1)
    s = iseven(M) ? -one(ρ) : one(ρ)
    L = (a / (ζ^2 + β2))^M                    # aᴹ·(ζ²+β²)^{-M}, |base|≤1 ⇒ underflows cleanly
    return _residue_Hm(iβ, ζ, s * SD1, -s * dζ * SD2, s * dζ^2 * SD3, L, σ)
end

# Non-integer M: residue fails at the branch point.
# H₀ is the Mace–Hellberg kappa-Z, a single Gauss ₂F₁ (Euler integral G&R 3.259.3)
function _kappa_Hm(ζ, β2, M, σ = 1)
    H₀ = _kappa_H0(ζ, β2, M, σ)
    N₀ = sqrt(π) * exp(loggamma(M - 0.5) - loggamma(M)) * β2^((1 - 2M) / 2)
    H₁ = N₀ + ζ * H₀
    H₂ = ζ * H₁
    return H₀, H₁, H₂
end

_kappa_Hm_scaled(ζ, a, b, M, σ = 1) = a^M .* _kappa_Hm(ζ, a * b, M, σ)

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
function para_moments(p::Kappa, Δ, kz)
    κ, a = p.kappa, p.a
    if iszero(kz)
        # no u-pole: (1, 0, ⟨u²⟩=a/(2κ−1), 0, ∫uf′=−1)/Δ
        invΔ = 1 / Δ
        return (invΔ, zero(invΔ), a / (2κ - 1) * invΔ, zero(invΔ), -invΔ)
    end
    ζ = Δ / kz
    σ = sign(kz)
    aH0, aH1, aH2 = _kappa_Hm_scaled(ζ, a, one(a), κ + 1, σ)
    _, aG1, aG2 = _kappa_Hm_scaled(ζ, a, one(a), κ + 2, σ)
    Cn = exp(loggamma(κ + 1) - loggamma(κ + 0.5)) / sqrt(π * a)     # 1-D norm
    pf = -Cn / kz
    tf = 2 * Cn * (κ + 1) / (a * kz)
    return (pf * aH0, pf * aH1, pf * aH2, tf * aG1, tf * aG2)
end

# uses the shared fused single-pass loop
function _separable_harmonics(para, p::Kappa, args...; kw...)
    κ, a = p.kappa, p.a
    C = κ / (π * a)                          # 2-D normalization

    f = v -> C * (1 + v^2 / a)^(-(κ + 1))
    fdf = v -> (D = 1 + v^2 / a; (C * D^(-(κ + 1)), -2C * (κ + 1) * v / a * D^(-(κ + 2))))
    fperp = AnalyticFactor{typeof(a)}(f, fdf)
    return _separable_harmonics_sum_first(para, fperp, args...; kw...)
end

# Slice moments Gₘ=∫uᵐ(b+u²/a∥)^{-M}/(u-ζ)du = a∥^M·H_m(ζ, a∥b, M), M=κ+2.
@inline _kappa_Gm(ζ, a_para, b, κ, σ = 1) = _kappa_Hm_scaled(ζ, a_para, b, κ + 2, σ)

# kz=0 slice: pole-free plain moments a∥^M·Sₘ=∫uᵐ(b+u²/a∥)^{-M}du (S₁=0 by parity), M=κ+2.
@inline function _kappa_Gm0(a_para, b, κ)
    M = κ + 2
    lgM = loggamma(M)
    S0 = sqrt(π) * exp(loggamma(M - 0.5) - lgM) * sqrt(a_para) * b^(0.5 - M)
    S2 = sqrt(π) / 2 * exp(loggamma(M - 1.5) - lgM) * a_para^1.5 * b^(1.5 - M)
    return S0, S2
end
