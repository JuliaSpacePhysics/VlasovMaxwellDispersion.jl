"""
    SeparableVDF(fpar, fperp; parlower, parupper, perpupper, dfpar=nothing, dfperp=nothing)
    SeparableVDF(fpar; lower, upper, df=nothing)

Arbitrary **separable analytic** VDF `f(v∥,v⊥) = f∥(v∥)·f⊥(v⊥)` for the **full
magnetized EM** susceptibility at oblique propagation (`k⊥≠0`). Both factors must
be evaluable at complex argument (continued onto the Landau contour).

The one-argument form is a reduced parallel distribution for the field-aligned
electrostatic path (`k⊥=0`): Landau damping / two-stream / bump-on-tail for any
analytic `f∥`.

Parallel moments close via the generic `hilbert` primitive; perp Bessel moments
by adaptive quadrature. The 3×3 tensor algebra is shared with the bi-Maxwellian
path (validated: `SeparableVDF(Gaussian,Gaussian) == Maxwellian`).
"""
struct SeparableVDF{Fp,Dp,Fq,Dq,T,Q} <: AbstractVDF
    fpar::Fp
    dfpar::Dp
    fperp::Fq
    dfperp::Dq
    parlo::T
    parhi::T
    perphi::Q
end

@inline reduced(d::SeparableVDF) = isnothing(d.fperp)

function SeparableVDF(fpar; lower, upper, df=nothing, normalize=true)
    lo, up = promote(float(lower), float(upper))
    n = normalize ? QuadGK.quadgk(fpar, lo, up; rtol=1.0e-10)[1] : one(lo)
    fp = u -> fpar(u) / n
    dfp = isnothing(df) ? (u -> _dwrt(fp, u)) : (u -> df(u) / n)
    return SeparableVDF(fp, dfp, nothing, nothing, lo, up, nothing)
end

function SeparableVDF(
    fpar, fperp; parlower, parupper, perpupper,
    dfpar=nothing, dfperp=nothing, normalize=true
)
    plo, phi = promote(float(parlower), float(parupper))
    qhi = oftype(phi, perpupper)
    np = normalize ? QuadGK.quadgk(fpar, plo, phi; rtol=1.0e-10)[1] : one(plo)
    nq = normalize ? 2π * QuadGK.quadgk(v -> fperp(v) * v, zero(qhi), qhi; rtol=1.0e-10)[1] : one(plo)
    fp = u -> fpar(u) / np
    fq = v -> fperp(v) / nq
    dfp = isnothing(dfpar) ? (u -> _dwrt(fp, u)) : (u -> dfpar(u) / np)
    dfq = isnothing(dfperp) ? (v -> _dwrt(fq, v)) : (v -> dfperp(v) / nq)
    return SeparableVDF(fp, dfp, fq, dfq, plo, phi, qhi)
end

# --- Arbitrary separable analytic f, full magnetized EM (oblique k⊥≠0) --------
# Same harmonic algebra as the bi-Maxwellian, but moments are computed generically:
# parallel z*F/z*T via the analytic `hilbert`, perp Bessel moments by quadrature.
function contribution(d::SeparableVDF, s::Species, ω, k::Wavenumber; kwargs...)
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    if iszero(kperp)
        reduced(d) && return _reduced_electrostatic_contribution(d, s, ω, k)
        throw(
            ArgumentError(
                "SeparableVDF: kperp=0 full EM tensor is not implemented; use reduced SeparableVDF(f; lower, upper) with electrostatic_det"
            )
        )
    end
    reduced(d) && throw(
        ArgumentError(
            "SeparableVDF: reduced one-argument form only supports field-aligned electrostatic kperp=0"
        )
    )
    ω = complex(float(ω))
    a = kperp / Ω                                   # k⊥/Ω, uniform Bessel arg coeff
    v⊥²_mean = 2π * QuadGK.quadgk(v -> d.fperp(v) * v^3, zero(d.perphi), d.perphi; rtol=1.0e-8)[1]
    nmax = nmax_bessel(a^2 * v⊥²_mean / 2)          # harmonic cap from the perp scale
    f = n -> _separable_harmonic(n, d, ω, Ω, kz, kperp, a)
    χ = converge(f, 1, 1.0e-7; nmax)
    return SMatrix{3,3,ComplexF64}((s.Pi2 / ω^2) * χ)
end

# χ_zz = -(Π²/k∥²) ∫ f∥′(u)/(u − ω/k∥) du
# Returns diag(0,0,χ_zz)
function _reduced_electrostatic_contribution(d::SeparableVDF, s::Species, ω, k::Wavenumber)
    kz = para(k)
    ω = complex(float(ω))
    χzz = -(s.Pi2 / kz^2) * hilbert(d.dfpar, ω / kz; lower=d.parlo, upper=d.parhi)
    z = zero(χzz)
    return @SMatrix ComplexF64[z z z; z z z; z z χzz]
end

@inline _besselj_prime(m, x) = (besselj(m - 1, x) - besselj(m + 1, x)) / 2

function _separable_harmonic(n, d::SeparableVDF, ω, Ω, kz, kperp, a)
    ζ = (ω - n * Ω) / kz
    L, U = d.parlo, d.parhi
    # Parallel: Landau–Hilbert for [f∥, u·f∥, u²·f∥, f∥′, u·f∥′]; the −1/kz folds the resonance kz.
    gpar(u) = (fp=d.fpar(u); dp=d.dfpar(u); SVector(fp, u * fp, u^2 * fp, dp, u * dp))
    z = (-1 / kz) .* hilbert(gpar, ζ; lower=L, upper=U)
    z0F, z1F, z2F, z0T, z1T = z[1], z[2], z[3], z[4], z[5]
    # Perp: Bessel-moment quadrature over [0, perphi] (∂Jn ≡ Jn′ wrt argument)
    Q = d.perphi
    function perp6(v)
        Jn = besselj(n, a * v);
        Jn′ = _besselj_prime(n, a * v)
        fq = d.fperp(v);
        dfq = d.dfperp(v)
        SVector(Jn^2 * fq * v, Jn^2 * dfq, Jn * Jn′ * dfq * v,
            Jn * Jn′ * fq * v^2, Jn′^2 * dfq * v^2, Jn′^2 * fq * v^3)
    end
    P = 2π .* QuadGK.quadgk(perp6, zero(Q), Q; rtol=1.0e-8)[1]
    JF, J∂F, JdJ∂F, JdJF, ∂J²∂F, ∂J²F = P[1], P[2], P[3], P[4], P[5], P[6]
    p = (; JF, J∂F, JdJF, JdJ∂F, ∂J²F, ∂J²∂F)
    return _chi_mblock((z0F, z1F, z2F, z0T, z1T), p, ω, kz, kperp, n / a)
end
