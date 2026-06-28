"""
    SeparableVDF(fperp, fpar; parlower, parupper, perpupper, dfpar=nothing, dfperp=nothing)
    SeparableVDF(fpar; lower, upper, df=nothing)

Arbitrary **separable analytic** VDF `f(p⊥,p∥) = f⊥(p⊥)·f∥(p∥)` for the **full
magnetized EM** susceptibility at oblique propagation (`k⊥≠0`). Both must be evaluable at complex arguments (continued onto the Landau contour).

The one-argument form is a reduced parallel distribution for the field-aligned
electrostatic path (`k⊥=0`): Landau damping / two-stream / bump-on-tail.
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

# Re-derive the oblique tensor from a closed-form Separable by quadrature over its callable factors.
SeparableVDF(d::Separable; kwargs...) = SeparableVDF(d.fperp, d.fpar; kwargs...)

function SeparableVDF(
    fperp, fpar; parlower, parupper, perpupper,
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
function contribution(d::SeparableVDF, s, ω, k; kwargs...)
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    reduced(d) && (iszero(kperp) ? (return _reduced_electrostatic_contribution(d, s, ω, k)) :
        throw(ArgumentError("SeparableVDF: reduced one-argument form only supports field-aligned electrostatic kperp=0")))
    ω = complex(float(ω))
    β = kperp / Ω                                   # k⊥/Ω, uniform Bessel arg coeff
    v⊥²_mean = 2π * QuadGK.quadgk(v -> d.fperp(v) * v^3, zero(d.perphi), d.perphi; rtol=1.0e-8)[1]
    nmax = nmax_bessel(β^2 * v⊥²_mean / 2)          # harmonic cap from the perp scale
    f = n -> _separable_harmonic(n, d, ω, Ω, kz, β)
    χ = converge(f, 1, 1.0e-7; nmax)
    return (s.Pi2 / ω^2) * χ
end

# χ_zz = -(Π²/k∥²) ∫ f∥′(u)/(u − ω/k∥) du
# Returns diag(0,0,χ_zz)
function _reduced_electrostatic_contribution(d::SeparableVDF, s, ω, k)
    kz = para(k)
    ω = complex(float(ω))
    χzz = -(s.Pi2 / kz^2) * hilbert(d.dfpar, ω / kz; lower=d.parlo, upper=d.parhi)
    z = zero(χzz)
    return @SMatrix ComplexF64[z z z; z z z; z z χzz]
end

@inline _besselj_prime(m, x) = (besselj(m - 1, x) - besselj(m + 1, x)) / 2

function _separable_harmonic(n, d::SeparableVDF, ω, Ω, kz, a)
    ζ = (ω - n * Ω) / kz
    L, U = d.parlo, d.parhi
    # Parallel: Landau–Hilbert for [f∥, u·f∥, u²·f∥, f∥′, u·f∥′]; the −1/kz folds the resonance kz.
    gpar(u) = (fp=d.fpar(u); dp=d.dfpar(u); SVector(fp, u * fp, u^2 * fp, dp, u * dp))
    M = (-1 / kz) .* hilbert(gpar, ζ; lower=L, upper=U)
    # Perp moments by quadrature over [0, perphi]: the 6 distinct entries of the
    # symmetric kernel bvec⊗bvec, bvec=(v⊥Rn, v⊥Jn′, Jn) with the ring kernel
    # Rn=(n/z)Jn=½(J_{n−1}+J_{n+1}) (regular at k⊥=0), in each of the f⊥′ (∂F) and
    # f⊥ (F) slices. F carries one extra v⊥.
    Q = d.perphi
    function perptri(v)
        z = a * v
        Jm, Jp1 = besselj(n - 1, z), besselj(n + 1, z)
        bvec = SVector(v * (Jm + Jp1) / 2, v * (Jm - Jp1) / 2, besselj(n, z))
        k11, k12, k13, k22, k23, k33 =
            bvec[1]^2, bvec[1] * bvec[2], bvec[1] * bvec[3], bvec[2]^2, bvec[2] * bvec[3], bvec[3]^2
        dfq, vfq = d.dfperp(v), v * d.fperp(v)
        SVector(dfq * k11, dfq * k12, dfq * k13, dfq * k22, dfq * k23, dfq * k33,
            vfq * k11, vfq * k12, vfq * k13, vfq * k22, vfq * k23, vfq * k33)
    end
    P = 2π .* QuadGK.quadgk(perptri, zero(Q), Q; rtol=1.0e-8)[1]
    P∂ = _symmat(P[1], P[2], P[3], P[4], P[5], P[6])
    PF = _symmat(P[7], P[8], P[9], P[10], P[11], P[12])
    return _chi_mblock(M, P∂, PF, ω, kz, n * Ω)
end
