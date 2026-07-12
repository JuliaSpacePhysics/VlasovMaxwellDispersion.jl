"""
    ReducedVDF(fpar; para=(lo,hi), df=nothing)

Reduced 1-D parallel distribution for the field-aligned electrostatic path (`k⊥=0`):
Landau damping / two-stream / bump-on-tail. `χ_zz = −(Π²/k∥²)∫ f∥′(u)/(u − ω/k∥) du`,
with `f∥` evaluable at complex argument.
"""
struct ReducedVDF{D, T, N} <: AbstractVDF
    df::D
    para::T
    n::N
end

function ReducedVDF(fpar; para, df = nothing, normalize = true)
    lo, hi = promote(float(para[1]), float(para[2]))
    n = normalize ? QuadGK.quadgk(fpar, lo, hi; rtol = 1.0e-10)[1] : one(lo)
    dfp = @something df (u -> _dwrt(fpar, u))
    return ReducedVDF(erase_f1(dfp, hi), (lo, hi), n)
end

function contribution(d::ReducedVDF, s, ω, k; rtol = 1.0e-9, closure = HarmonicSum(), kw...)
    iszero(perp(k)) ||
        throw(ArgumentError("ReducedVDF (1-D parallel) only supports field-aligned electrostatic kperp=0"))
    kz = para(k)
    χzz = -(s.Pi2 / kz^2) * plan_landau(d.para, ω / kz, sign(kz))(d.df; rtol, kw...)
    return _ee33(χzz) / d.n
end
