# Duck-typed contract (no factor supertypes). 
# A *parallel* factor defines
#   par_moments(par, ω, kz, nΩ)   -> M = (MF0,MF1,MF2,MT0,MT1)
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
    fpar::P
end

⊗(perp, par) = Separable(perp, par)
perp_setup(perp, β) = perp

(d::Separable)(q, u) = d.fperp(q) * d.fpar(u)

function contribution(d::Separable, s, ω, k; rtol = 1.0e-8, kwargs...)
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    ω = complex(float(ω))
    β = kperp / Ω
    return (s.Pi2 / ω^2) * _separable_harmonics(d.fpar, perp_setup(d.fperp, β), β, ω, Ω, kz, rtol)
end

# Function barrier: `prepared` type is value-dependent
function _separable_harmonics(par, perp, β, ω, Ω, kz, rtol)
    return converge(rtol; nmax = nmax_harm(perp, β)) do n
        nΩ = n * Ω
        M = par_moments(par, ω, kz, nΩ)
        P∂, PF = perp_moments(perp, n, β)
        return _chi_mblock(M, P∂, PF, ω, kz, nΩ)
    end
end
