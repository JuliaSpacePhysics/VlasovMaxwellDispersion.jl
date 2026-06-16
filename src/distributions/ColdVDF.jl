"""
    ColdVDF()

`T → 0` limit; handled via the closed cold (Stix S,D,P) form for accuracy.
"""
struct ColdVDF <: AbstractVDF end

Continuation(::ColdVDF) = Analytic()

@inline thermal_par(::ColdVDF) = false   # exact zero, type-stable additive identity
@inline thermal_perp(::ColdVDF) = false
@inline drift(::ColdVDF) = false

#   χ = ε_cold - I:
#   χ_xx = χ_yy = -Π²/(ω²-Ω²)
#   χ_xy = -χ_yx = -i Ω Π²/(ω(ω²-Ω²))   (sign tracks signed Ω ⇒ charge)
#   χ_zz = -Π²/ω²
function contribution(::ColdVDF, s, ω, k; kwargs...)
    Π2 = s.Pi2
    Ω = s.Omega
    ω = complex(float(ω))
    den = ω^2 - Ω^2
    S = -Π2 / den
    D = -im * Ω * Π2 / (ω * den)
    P = -Π2 / ω^2
    z = zero(S)
    return @SMatrix ComplexF64[
        S   D   z;
        -D  S   z;
        z   z   P
    ]
end
