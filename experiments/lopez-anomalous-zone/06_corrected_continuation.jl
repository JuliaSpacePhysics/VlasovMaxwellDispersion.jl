# The FIX for López's continuation: their θ term freezes the resonant support
# at γ1j(Re z); the analytic continuation instead continues the support
# endpoint into complex γ,
#   γ1e(z) = ( z t + √(t²+1−z²)) / (1−z²)      (electron)
#   γ1p(z) = (−z t + √(t²+1−z²)) / (1−z²)      (positron)
# and adds 2πi·σ/2 times the tail integral along a complex path from γ1j(z)
# to +∞ (integrand e^{−μγ}·w(γ) is entire, path is free). On the real axis
# this equals López exactly; off it, it is holomorphic where theirs is not.
# Not restricted to |Re z| < 1: the branch functions are analytic on any
# lower-half-plane path avoiding z = ±1 (the √ principal cut needs real
# z² > 1+t², the endpoint denominators pole at z = ±1), so the same formula
# followed continuously through the gap 1 < Re z < √(1+t²) is the unique
# continuation of the subluminal germ — see 09_superluminal_continuation.jl.
include("lopez.jl")
using Printf

γ1e_c(z, t) = (z * t + sqrt(t^2 + 1 - z^2)) / (1 - z^2)
γ1p_c(z, t) = (-z * t + sqrt(t^2 + 1 - z^2)) / (1 - z^2)

# ∫_{γ0}^{∞} e^{-μγ} w(γ) dγ, complex γ0: leg γ0→(real(γ0)+ i0-free landing)…
# integrand entire ⇒ straight leg γ0→a then real ray a→∞, a = max(real(γ0),1)+0.5
function tail(w, γ0, μ)
    a = max(real(γ0), 1.0) + 0.5
    leg1, _ = quadgk(s -> exp(-μ * (γ0 + s * (a - γ0))) * w(γ0 + s * (a - γ0)) * (a - γ0), 0.0, 1.0; rtol = 1.0e-10)
    leg2, _ = quadgk(γ -> exp(-μ * γ) * w(γ), a, Inf; rtol = 1.0e-10)
    return leg1 + leg2
end

# corrected Λ_L: θ-less principal part + analytic continuation term (Im z ≤ 0)
function ΛL_corr(x, y, μ)
    z = x / y
    t = 1 / y
    base = ΛL(
        x, y, μ;
        Jefun = (γ, zz, tt) -> (R = real(zz); I = imag(zz); s1 = S1e(γ, tt); s2 = S2e(γ, tt);
        complex(0.5 * log(((R - s2)^2 + I^2) / ((R + s1)^2 + I^2)), atan((s2 - R) / I) + atan((s1 + R) / I))),
        Jpfun = (γ, zz, tt) -> (R = real(zz); I = imag(zz); s1 = S1p(γ, tt); s2 = S2p(γ, tt);
        complex(0.5 * log(((R - s2)^2 + I^2) / ((R + s1)^2 + I^2)), atan((s2 - R) / I) + atan((s1 + R) / I)))
    )
    imag(z) > 0 && return base
    σ = imag(z) == 0 ? 1.0 : 2.0
    K2 = besselk(2, μ)
    pref = (μ^2 / (4 * K2)) / (x * y^3)
    we = γ -> (y^2 - x^2) * γ^2 - 2 * x * γ - (1 + y^2)
    wp = γ -> (y^2 - x^2) * γ^2 + 2 * x * γ - (1 + y^2)
    cont = im * π * σ * (tail(we, γ1e_c(z, t), μ) + tail(wp, γ1p_c(z, t), μ))
    return base + pref * cont
end

function dzbar(f, z; h = 1.0e-4)
    fx = (f(z + h) - f(z - h)) / (2h)
    fy = (f(z + im * h) - f(z - im * h)) / (2h)
    return (fx + im * fy) / 2, (fx - im * fy) / 2
end

println("== 1) Corrected = López on/near the real axis (continuity preserved) ==")
for x in (0.05, 0.124, 0.3)
    a, b = ΛL_corr(complex(x, -1.0e-8), 0.5, 2.0), ΛL(complex(x, -1.0e-8), 0.5, 2.0)
    @printf("  x=%.3f-1e-8im: |corr-López|/|López| = %.1e\n", x, abs(a - b) / abs(b))
end

println("\n== 2) Holomorphy |∂z̄|/|∂z| in LHP: corrected vs López-with-θ ==")
for (lbl, z) in [
        ("LHP  0.124-0.06i", 0.124 - 0.06im), ("desc-root 0.124-0.136i", 0.124 - 0.136im),
        ("LHP  0.124-0.25i", 0.124 - 0.25im),
    ]
    zbC, zC = dzbar(ω -> ΛL_corr(ω, 0.5, 2.0), z)
    zbL, zL = dzbar(ω -> ΛL(ω, 0.5, 2.0), z)
    @printf("  %-24s corrected=%.2e   López=%.2e\n", lbl, abs(zbC) / abs(zC), abs(zbL) / abs(zL))
end

println("\n== 3) Roots of the corrected continuation (Muller) ==")
println("   VMD reference roots: k=0.5: 0.16297-0.11048im ; k=0.6: 0.20540-0.18512im")
for (k, seeds) in ((0.5, (0.124 - 0.136im, 0.163 - 0.110im)), (0.6, (0.004 - 0.234im, 0.205 - 0.185im)))
    f = ω -> ΛL_corr(ω, k, 2.0)
    for s in seeds
        r = muller(f, s - 1.0e-3, s, s + 1.0e-3im)
        @printf(
            "  k=%.1f seed %+.3f%+.3fim -> root %+.5f%+.5fim  |Λ|=%.1e\n",
            k, real(s), imag(s), real(r), imag(r), abs(f(r))
        )
    end
end
println("\nExpected: both seeds converge to the VMD root; no root remains at the López descent location.")
