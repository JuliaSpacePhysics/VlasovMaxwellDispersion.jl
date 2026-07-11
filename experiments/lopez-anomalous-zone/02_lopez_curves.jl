include("lopez.jl")
using Printf

# Trace López Λ_L root in k from ALPS seed; does ωr descend (their fig) or rise (VMD)?
function trace(μ, ks, seed)
    out = ComplexF64[]
    s = seed
    for k in ks
        f = ω -> ΛL(ω, k, μ)
        r = muller(f, s-1e-3, s, s+1e-3im)
        push!(out, r)
        s = r
    end
    out
end

println("== μ=2: trace López Λ_L propagating root from ALPS seed ==")
μ=2.0
ks = 0.1:0.05:0.85
seed = complex(0.039185,-2.5e-5)
rts = trace(μ, ks, seed)
for (k,r) in zip(ks,rts)
    @printf("  k=%.2f  ω=%+.5f%+.5eim  |Λ|=%.1e\n", k, real(r), imag(r), abs(ΛL(r,k,μ)))
end

println("\n== μ=2: seed directly on digitized descent (ωr≈0.08→0.02, γ≈-0.15→-0.3) ==")
descseed = [(0.50,complex(0.066,-0.15)),(0.55,complex(0.066,-0.22)),
            (0.60,complex(0.04,-0.24)),(0.65,complex(0.022,-0.27))]
for (k,s) in descseed
    f = ω -> ΛL(ω, k, μ)
    r = muller(f, s-1e-2, s, s+1e-2im)
    @printf("  k=%.2f seed=%+.3f%+.3fim -> ω=%+.5f%+.5fim |Λ|=%.1e\n",
        k, real(s),imag(s), real(r),imag(r), abs(f(r)))
end

println("\n== μ=10: trace from k=0.3 seed upward through the peak+descent ==")
μ=10.0
ks10 = 0.3:0.1:2.6
seed10 = complex(0.155,-1e-3)
r10 = trace(μ, ks10, seed10)
for (k,r) in zip(ks10,r10)
    @printf("  k=%.2f  ω=%+.5f%+.5fim  |Λ|=%.1e\n", k, real(r), imag(r), abs(ΛL(r,k,μ)))
end

println("\n== μ=10: seed directly on digitized descent ==")
d10 = [(1.90,complex(0.42,-0.72)),(2.15,complex(0.275,-0.98)),(2.45,complex(0.077,-1.19))]
for (k,s) in d10
    f = ω -> ΛL(ω, k, μ)
    r = muller(f, s-1e-2, s, s+1e-2im)
    @printf("  k=%.2f seed=%+.3f%+.3fim -> ω=%+.5f%+.5fim |Λ|=%.1e\n",
        k, real(s),imag(s), real(r),imag(r), abs(f(r)))
end
