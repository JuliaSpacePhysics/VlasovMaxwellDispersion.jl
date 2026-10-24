# Payoff check for theory.md §3: does the gyroharmonic ladder move as eq. (3) says, and
# does that move the current threshold? Topside auroral O⁺/e⁻, σ kept inside the §5 scope.

include(joinpath(@__DIR__, "shear.jl"))
using Printf

const MI_ME = 16 * 1836.15
const EPS = 6.321952600232268e-6     # w(O⁺, 0.3 eV)/c
const PI2 = 11565.761430014622       # (ω_pO⁺/Ω_O⁺)², n = 3×10³ cm⁻³, B = 2.8×10⁻⁵ T
const KRHO, RATIO = 1.5, 0.05        # peak of the unsheared current-driven EIC band

kvec(krho = KRHO, r = RATIO) = ShearK(0.0, krho / EPS, r * krho / EPS)
plasma(u) = (NormalizedSpecies(1.0, PI2, Maxwellian(vth = EPS)),
    NormalizedSpecies(-MI_ME, PI2 * MI_ME,
        Maxwellian(vth = EPS * sqrt(MI_ME), vd = u * EPS * sqrt(MI_ME))))

function muller(f, x0, x1, x2; tol = 1e-11, maxit = 60)
    f0, f1, f2 = f(x0), f(x1), f(x2)
    for _ in 1:maxit
        q = (x2 - x1) / (x1 - x0)
        A = q * f2 - q * (1 + q) * f1 + q^2 * f0
        B = (2q + 1) * f2 - (1 + q)^2 * f1 + q^2 * f0
        C = (1 + q) * f2
        d = sqrt(B^2 - 4A * C)
        den = abs(B + d) > abs(B - d) ? B + d : B - d
        iszero(den) && break
        x3 = x2 - (x2 - x1) * 2C / den
        x0, x1, x2, f0, f1, f2 = x1, x2, x3, f1, f2, f(x3)
        abs(x2 - x1) < tol * max(1, abs(x2)) && break
    end
    return x2
end

eic_root(u, S; seed = 1.185 + 0.02im, k = kvec()) =
    muller(ω -> shear_epsilon(plasma(u), ω, k, S), seed * 0.98, seed * 1.02, seed)

"Bisect the electron drift whose EIC growth vanishes."
function threshold_u(S; lo = 0.0, hi = 0.6, iters = 22)
    for _ in 1:iters
        mid = (lo + hi) / 2
        imag(eic_root(mid, S)) > 0 ? (hi = mid) : (lo = mid)
    end
    return hi
end

println("EIC branch vs shear, O⁺/e⁻, k⊥ρ = $KRHO, k∥/k⊥ = $RATIO, electron drift u = 0.2 v_e\n")
@printf("%7s  %8s  %9s  %9s  %10s  %8s\n", "σ", "√(1+σ)", "ω_r/Ω", "γ/Ω", "u_c/v_e", "u_c/u_c(0)")
u0 = threshold_u(0.0)
for σ in (-0.10, -0.05, -0.02, 0.0, 0.02, 0.05, 0.10)
    ω = eic_root(0.2, σ)
    uc = threshold_u(σ)
    @printf("%+7.3f  %8.4f  %9.4f  %+9.5f  %10.4f  %8.3f\n",
        σ, sqrt(1 + σ), real(ω), imag(ω), uc, uc / u0)
end
