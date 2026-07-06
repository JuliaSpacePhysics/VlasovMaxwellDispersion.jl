"""
    ReducedVDF(fpar; para=(lo,hi), df=nothing)

Reduced 1-D parallel distribution for the field-aligned electrostatic path (`k⊥=0`):
Landau damping / two-stream / bump-on-tail. `χ_zz = −(Π²/k∥²)∫ f∥′(u)/(u − ω/k∥) du`,
with `f∥` evaluable at complex argument.
"""
struct ReducedVDF{D, T} <: AbstractVDF
    df::D       # normalized f∥′
    para::T
end

function ReducedVDF(fpar; para, df = nothing, normalize = true)
    lo, hi = promote(float(para[1]), float(para[2]))
    n = normalize ? QuadGK.quadgk(fpar, lo, hi; rtol = 1.0e-10)[1] : one(lo)
    fp = u -> fpar(u) / n
    dfp = isnothing(df) ? (u -> _dwrt(fp, u)) : (u -> df(u) / n)
    return ReducedVDF(dfp, (lo, hi))
end

function contribution(d::ReducedVDF, s, ω, k; kwargs...)
    iszero(perp(k)) ||
        throw(ArgumentError("ReducedVDF (1-D parallel) only supports field-aligned electrostatic kperp=0"))
    kz = para(k)
    χzz = -(s.Pi2 / kz^2) * hilbert(d.df, ω / kz, d.para...; σ = sign(kz))
    z = zero(χzz)
    return @SMatrix ComplexF64[z z z; z z z; z z χzz]
end
