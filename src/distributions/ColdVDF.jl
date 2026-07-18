"""
    ColdVDF()

`T → 0` limit; handled via the closed cold (Stix S,D,P) form for accuracy.
"""
struct ColdVDF <: AbstractVDF end
parallel_even(::ColdVDF) = true

# Gyrotropic tensor [S D 0; -D S 0; 0 0 P]
_gyrotropic(S, D, P) = (z = zero(complex(S)); @SMatrix ComplexF64[S D z; -D S z; z z P])

#   χ = ε_cold - I:  χ_xx=χ_yy = -Π²/(ω²-Ω²),  χ_xy = -iΩΠ²/(ω(ω²-Ω²)),  χ_zz = -Π²/ω²
function contribution(::ColdVDF, s, ω, k; kwargs...)
    Π2, Ω, den = s.Pi2, s.Omega, ω^2 - s.Omega^2
    return _gyrotropic(-Π2 / den, -im * Ω * Π2 / (ω * den), -Π2 / ω^2)
end

# ω̃²·χ with the 1/ω² (P) and 1/ω (D) poles cancelled analytically
function scaled_contribution(::ColdVDF, s, ω, k; kwargs...)
    Π2, Ω = s.Pi2, s.Omega
    den = ω^2 - Ω^2
    return _gyrotropic(-Π2 * ω^2 / den, -im * Ω * Π2 * ω / den, -Π2)
end
