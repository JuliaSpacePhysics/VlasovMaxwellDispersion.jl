include("lopez.jl")
using Printf
ΛL_noθ(x,y,μ) = ΛL(x,y,μ;
  Jefun=(γ,z,t)->(R=real(z);I=imag(z);s1=S1e(γ,t);s2=S2e(γ,t);
    complex(0.5*log(((R-s2)^2+I^2)/((R+s1)^2+I^2)), atan((s2-R)/I)+atan((s1+R)/I))),
  Jpfun=(γ,z,t)->(R=real(z);I=imag(z);s1=S1p(γ,t);s2=S2p(γ,t);
    complex(0.5*log(((R-s2)^2+I^2)/((R+s1)^2+I^2)), atan((s2-R)/I)+atan((s1+R)/I))))
function dzbar(f,z;h=1e-4)
    fx=(f(z+h)-f(z-h))/(2h); fy=(f(z+im*h)-f(z-im*h))/(2h); (fx+im*fy)/2,(fx-im*fy)/2
end
println("== Mechanism: holomorphy of Lopez WITHOUT θ term (θ isolated as sole non-analyticity) ==")
for (lbl,z) in [("UHP 0.124+0.30i",0.124+0.30im),("LHP 0.124-0.06i",0.124-0.06im),
                ("desc-root 0.124-0.136i",0.124-0.136im),("LHP 0.124-0.25i",0.124-0.25im)]
    zb,zz=dzbar(ω->ΛL_noθ(ω,0.5,2.0),z)
    @printf("  %-24s |∂z̄|/|∂z| (no-θ)=%.2e\n",lbl,abs(zb)/abs(zz))
end
f = ω -> ΛL(ω, 0.5, 2.0)
r = muller(f, 0.123-0.135im, 0.124-0.136im, 0.125-0.137im)
@printf("\nLopez descent root (Muller-refined): ω=%.6f%+.6fim, |ΛL| = %.2e — genuine zero of the\nnon-holomorphic fn (steep: |ΛL| at 5-digit rounding is already %.2e)\n",
    real(r), imag(r), abs(f(r)), abs(ΛL(0.12415-0.13610im,0.5,2.0)))
