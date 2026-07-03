"""
    BiKappa(; vth_para, vth_perp=vth_para, kappa)

Coupled anisotropic bi-kappa `f₀ ∝ (1 + p∥²/a∥ + p⊥²/a⊥)^{-(κ+1)}`,
`a_{∥,⊥} = (κ-3/2)·vth_{∥,⊥}²` (so `⟨p⊥²⟩=vth_perp²`), `κ>3/2`.
"""
struct BiKappa{K, T} <: AbstractVDF
    kappa::K
    a_para::T
    a_perp::T
end

function BiKappa(; vth_para, vth_perp = vth_para, kappa)
    κ = float(kappa)
    kappa > 1.5 || throw(ArgumentError("BiKappa needs κ > 3/2 (finite ⟨p²⟩)"))
    a_para, a_perp = promote((κ - 1.5) * vth_para^2, (κ - 1.5) * vth_perp^2)
    return BiKappa(kappa, a_para, a_perp)
end

(d::BiKappa)(q, u) = (1 + u^2 / d.a_para + q^2 / d.a_perp)^(-(d.kappa + 1))

# ∫d³p f₀=1 for the un-normalized (1+…)^{-(κ+1)}; closed-form moments carry C explicitly.
@inline normalization(d::BiKappa) =
    gamma(d.kappa + 1) / (π^1.5 * d.a_perp * sqrt(d.a_para) * gamma(d.kappa - 0.5))

function contribution(d::BiKappa, s, ω, k; rtol = 1.0e-8, norm = NORM, kwargs...)
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    ω = complex(float(ω))
    a = kperp / Ω
    p⊥2 = d.a_perp / (d.kappa - 1.5)
    ns = (-nmax_bessel(a^2 * p⊥2 / 2)):nmax_bessel(a^2 * p⊥2 / 2)
    vc = sqrt(p⊥2)
    C = normalization(d)
    X = QuadGK.quadgk(zero(vc), oftype(vc, Inf); rtol, norm) do v
        _harmonic_sum_perp(d, v, ns, C, ω, Ω, kz, a)
    end[1]
    return (s.Pi2 / ω^2) * _antisymmat(X)
end

# at fixed p⊥ the parallel slice is a 1-D kappa,
#     f(·,p⊥) ∝ (b + p∥²/a∥)^{-(κ+1)},  b = 1 + p⊥²/a⊥,
# whose Cauchy moments close (β²=a∥b in H_m).
# invk is folded into z, so `_In_block`'s c=1 supplies only the ∫d³p 2π factor.
function _harmonic_sum_perp(d::BiKappa, v, ns, C, ω, Ω, kz, a)
    κ, a_para, a_perp = d.kappa, d.a_para, d.a_perp
    b = 1 + v^2 / a_perp
    invk = -1 / kz
    # ∂⊥f = cF·D^{-M},  ∂∥f = cT·u·D^{-M},  M=κ+2, D=b+u²/a∥
    cF = invk * (-2 * C * (κ + 1) * v / a_perp)
    cT = invk * (-2 * C * (κ + 1) / a_para)
    return @no_escape begin
        b2s = @alloc(SVector{6, typeof(a * v)}, length(ns))
        _perp_Bessel_bilinears!(b2s, a, v)
        acc = zero(AType)
        @inbounds for (i, n) in enumerate(ns)
            ζ = (ω - n * Ω) / kz
            G0, G1, G2 = _kappa_Gm(ζ, a_para, b, κ)
            # F slice ← ∂⊥f (moments G0..G2), T slice ← ∂∥f (uᵐ·∂∥f → G_{m+1})
            z = (cF * G0, cF * G1, cF * G2, cT * G1, cT * G2)
            acc += _In_block(z, 1, b2s[i], v, ω, kz, n * Ω)
        end
        acc
    end
end
