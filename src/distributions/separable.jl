# Duck-typed contract (no factor supertypes).
# A *parallel* factor defines
#   para_moments(para, ω, kz, nΩ)   -> M = (MF0,MF1,MF2,MT0,MT1)
# A *perpendicular* factor defines
#   perp_setup(perp, β)           -> prepared   (default: itself; rings/tables override)
#   nmax_harm(prepared, β)        -> Int        (harmonic cap from the perp scale)
#   perp_moments(prepared, n, β)  -> (P∂, PF)

"""
    Separable(fperp, fpar)
    fperp ⊗ fpar

Closed-form separable VDF `f = fperp(p⊥)·fpar(p∥)` from 1D factors.
"""
struct Separable{Q, P} <: AbstractVDF
    fperp::Q
    fpara::P
end

⊗(perp, para) = Separable(perp, para)
perp_setup(perp, β) = perp

(d::Separable)(q, u) = d.fperp(q) * d.fpara(u)

function contribution(d::Separable, s, ω, k; rtol = 1.0e-8, kwargs...)
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    β = kperp / Ω
    X = _separable_harmonics(d.fpara, perp_setup(d.fperp, β), β, ω, Ω, kz; rtol)
    return _antisymmat((s.Pi2 / ω^2) * X)
end

function nmax_harm(p, β)
    p⊥²_mean = 2π * QuadGK.quadgk(v -> p(v) * v^3, p.lo, p.hi; rtol = 1.0e-8)[1]
    return nmax_bessel(β^2 * abs(p⊥²_mean) / 2)
end

# Function barrier: `prepared` type is value-dependent
function _separable_harmonics(para, perp, args...; kwargs...)
    return _separable_harmonics_sum_last(para, perp, args...; kwargs...)
end

function _separable_harmonics_sum_last(para, perp, β, ω, Ω, kz; rtol)
    return converge(; nmax = nmax_harm(perp, β), rtol) do n
        nΩ = n * Ω
        Δ = ω - nΩ
        M = para_moments(para, Δ, kz)
        P∂, PF = perp_moments(perp, n, β)
        return _chi_mblock(M, P∂, PF, ω, kz, nΩ)
    end
end

# Fused single-pass harmonic loop: parallel moments Mₙ are v-independent
function _separable_harmonics_sum_first(para, perp, β, ω, Ω, kz; rtol, norm = NORM)
    nmax = nmax_harm(perp, β)
    ns = -nmax:nmax
    Ms = _para_moments_all(para, ω, kz, Ω, ns; rtol)
    M = last(ns) + 1
    return @no_escape begin
        Jv = @alloc(typeof(β), M + 1)
        QuadGK.quadgk(perp.lo, perp.hi; rtol, norm) do v
            z = β * v
            fq, dfq = perp.fdf(v)
            vfq = v * fq
            besselj_ladder!(Jv, M, z)        # J_0..J_{nmax+1} in one recurrence, signed-indexed
            sum(enumerate(ns)) do (i, n)
                Jm, Jn, Jp = _jladder(Jv, n - 1), _jladder(Jv, n), _jladder(Jv, n + 1)
                # bvec=(v⊥Rn, v⊥Jn′, Jn), Rn=½(J_{n−1}+J_{n+1}); K=bvec⊗bvec shared by ∂F/F slices
                b1, b2, b3 = v * (Jm + Jp) / 2, v * (Jm - Jp) / 2, Jn
                K = _symmat(b1^2, b1 * b2, b1 * b3, b2^2, b2 * b3, b3^2)
                _chi_mblock(Ms[i], (2π * dfq) .* K, (2π * vfq) .* K, ω, kz, n * Ω)
            end
        end[1]
    end
end

_para_moments_all(p, ω, kz, Ω, ns; rtol = 1.0e-8) = map(n -> para_moments(p, ω - n * Ω, kz), ns)
