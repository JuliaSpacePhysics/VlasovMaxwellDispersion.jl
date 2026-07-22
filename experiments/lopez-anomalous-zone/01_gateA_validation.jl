include("lopez.jl")
using Printf, Random
Random.seed!(1)

println("== Gate A: closed-form J vs direct integral, Im z>0 (θ must be 0) ==")
maxerr = 0.0; maxth = 0.0
for _ in 1:20
    γ = 1 + 3*rand()
    z = complex(2*(rand()-0.5), 0.05 + rand())   # Im z>0
    t = 0.5 + 10*rand()
    # only meaningful where sqrt(1-1/γ²) real (γ>1) always ok
    de = Je(γ,z,t); dd = Je_direct(γ,z,t)
    dp = Jp(γ,z,t); ddp = Jp_direct(γ,z,t)
    global maxerr = max(maxerr, abs(de-dd), abs(dp-ddp))
    global maxth = max(maxth, abs(θe(γ,real(z),imag(z),t)), abs(θp(γ,real(z),imag(z),t)))
end
@printf("max |Je_closed - Je_direct| = %.2e ; max |θ| in UHP = %.2e\n", maxerr, maxth)

println("\n== Gate A2: Λ_L zero at ALPS root ω=0.039621-2.644e-6im, k=0.1, μ=2 ==")
μ=2.0; k=0.1
for lbl_ω in (("+ω", complex(3.9621e-2,-2.644e-6)), ("-conj(ω) mirror", -conj(complex(3.9621e-2,-2.644e-6))),
              ("VMD 0.03919", complex(0.03919,-2.6e-6)))
    lbl, ω = lbl_ω
    v = ΛL(ω, k, μ)
    @printf("  %-16s Λ_L=%+.4e%+.4eim  |Λ_L|=%.3e\n", lbl, real(v), imag(v), abs(v))
end
# refine the L-mode root near real axis with Muller
f = ω -> ΛL(ω, k, μ)
r = muller(f, 0.035+0.0im, 0.040+0.0im, 0.045-1e-4im)
@printf("  Muller root near 0.04: ω=%.6f%+.3eim  |Λ|=%.2e\n", real(r), imag(r), abs(f(r)))
