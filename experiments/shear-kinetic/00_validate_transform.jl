# Numerical arbiter for theory.md §3-§4. Two questions:
#   (a) is the orbit reduction exact?          → compare with kᵧu(X) switched off
#   (b) how big is the term it cannot absorb?  → error vs (σ, kᵧρ)
# Reference is the linear response integrated directly along exact sheared orbits in lab
# velocity coordinates: no invariants, no circularizing map, no Bessel expansion.

include(joinpath(@__DIR__, "shear.jl"))
using QuadGK, Printf
using LinearAlgebra: SymTridiagonal, eigen

"Gauss–Hermite nodes/weights for ∫f(x)e^{-x²}dx (Golub–Welsch)."
function gausshermite(n)
    E = eigen(SymTridiagonal(zeros(n), [sqrt(k / 2) for k in 1:(n - 1)]))
    return E.values, sqrt(π) .* abs2.(E.vectors[1, :])
end

"""ε_L for one Maxwellian species of width `w` under linear shear `S`, by direct orbit
integration (theory.md §4 eq. 12). `doppler=false` drops the kᵧu(X) term alone."""
function epsilon_bruteforce(Ω, Π2, w, k::ShearK, S, ω; n = 40, smax = 60.0, doppler = true)
    f = 1 + S / Ω
    Ωe = Ω * sqrt(f)
    t, W = gausshermite(n)
    acc = zero(complex(ω))
    for i in eachindex(t), j in eachindex(t)
        vx = w * t[i]
        vy = w * sqrt(f) * t[j]                 # lab weight is elliptical (theory.md §2)
        X = vy / (Ω * f)                        # guiding centre, u₀ = 0
        uX = S * X
        function integrand(s)
            c, sn = cos(Ωe * s), -sin(Ωe * s)   # orbit evaluated at τ = −s
            Δx = X * (1 - c) + (vx / Ωe) * sn
            Δy = (doppler ? -uX * s : zero(uX)) + (Ω * X / Ωe) * sn - (Ω * vx / Ωe^2) * (1 - c)
            vxt = X * Ωe * sn + vx * c
            vyt = Ω * X * c - (Ω * vx / Ωe) * sn        # vᵧ(τ) − u(X)
            return cis(k.kx * Δx + k.ky * Δy + ω * s) * exp(-k.kz^2 * w^2 * s^2 / 4) *
                   (-(2 / w^2) * (k.kx * vxt + k.ky * vyt) + im * k.kz^2 * s)
        end
        acc += W[i] * W[j] * first(quadgk(integrand, 0.0, 2.0, 6.0, 15.0, 30.0, smax;
            rtol = 1e-9, order = 15))
    end
    return 1 - im * (Π2 / abs2(k)) * acc / π
end

const W = 0.05
const OM = 1.2 + 0.3im                   # Im ω > 0 ⇒ absolutely convergent, no Landau contour
species() = (NormalizedSpecies(1.0, 1.0, Maxwellian(vth = W)),)
kvec(ky) = ShearK(10.0, ky, 3.0)

println("(a) orbit reduction — theory.md eq. (7), single Maxwellian ion, k=(10,20,3), kᵧρ=1.0\n")
@printf("%7s  %-26s  %-26s  %9s  %9s\n", "S", "brute force (full)", "mapping eq. (7)", "err full", "err no kᵧu")
for S in (0.0, 0.1, -0.15, 0.3, -0.4, 0.8)
    K = kvec(20.0)
    mp = shear_epsilon(species(), OM, K, S)
    bf = epsilon_bruteforce(1.0, 1.0, W, K, S, OM)
    nd = epsilon_bruteforce(1.0, 1.0, W, K, S, OM; doppler = false)
    @printf("%7.2f  %+11.8f%+11.8fim  %+11.8f%+11.8fim  %9.2e  %9.2e\n",
        S, real(bf), imag(bf), real(mp), imag(mp), abs(bf - mp) / abs(mp), abs(nd - mp) / abs(mp))
end

println("\n(b) the term eq. (7) cannot absorb: |ε_exact − ε_eq7|/|ε_eq7| vs (σ, kᵧρ)")
@printf("%8s", "σ \\ kᵧρ")
kys = (8.0, 16.0, 24.0, 32.0)
for ky in kys
    @printf("%10.1f", ky * W)
end
println()
for S in (0.05, 0.1, 0.2, 0.4)
    @printf("%8.2f", S)
    for ky in kys
        K = kvec(ky)
        e = abs(epsilon_bruteforce(1.0, 1.0, W, K, S, OM) - shear_epsilon(species(), OM, K, S)) /
            abs(shear_epsilon(species(), OM, K, S))
        @printf("%10.3f", e)
    end
    println()
end
