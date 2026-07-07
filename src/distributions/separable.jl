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
function _separable_harmonics(para, perp, β, ω, Ω, kz; rtol)
    return converge(; nmax = nmax_harm(perp, β), rtol) do n
        nΩ = n * Ω
        Δ = ω - nΩ
        M = para_moments(para, Δ, kz)
        P∂, PF = perp_moments(perp, n, β)
        return _chi_mblock(M, P∂, PF, ω, kz, nΩ)
    end
end
