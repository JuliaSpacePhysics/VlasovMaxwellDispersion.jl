include("lopez.jl")
using VlasovMaxwellDispersion
using VlasovMaxwellDispersion: DispersionFunction
using RationalFunctionApproximation: approximate, poles
using Printf

pair(vdf) = (NormalizedSpecies(1.0, 1.0, vdf), NormalizedSpecies(-1.0, 1.0, vdf))
plasma2 = pair(MaxwellJuttner(mu = 2.0)); plasma10 = pair(MaxwellJuttner(mu = 10.0))
vmd_defl(pl, kz) = DispersionFunction(pl, Wavenumber(0.0, kz); mode = :L, deflate = true)

# fixed Wirtinger d/dz̄, d/dz : y-derivative uses real spacing h
function dzbar(f, z; h=1e-4)
    fx = (f(z+h) - f(z-h))/(2h)
    fy = (f(z+im*h) - f(z-im*h))/(2h)
    (fx + im*fy)/2, (fx - im*fy)/2
end
println("== Holomorphy |∂z̄|/|∂z| (0 ⇒ analytic). Lopez ΛL (with θ) vs VMD L-mode, k=0.5 μ=2 ==")
for (lbl,z) in [("UHP  0.124+0.30i",0.124+0.30im),("axis 0.124+0.02i",0.124+0.02im),
                ("LHP  0.124-0.06i",0.124-0.06im),("desc-root 0.124-0.136i",0.124-0.136im),
                ("LHP  0.124-0.25i",0.124-0.25im)]
    zbL,zL = dzbar(ω->ΛL(ω,0.5,2.0), z)
    zbV,zV = dzbar(vmd_defl(plasma2,0.5), z)
    @printf("  %-24s  Lopez=%.2e   VMD=%.2e\n", lbl, abs(zbL)/abs(zL), abs(zbV)/abs(zV))
end

# AAA: report, for each candidate, the closest extrapolated pole
function aaa_poles_res(f; rebox, imbox, nre, nim)
    Z=[complex(x,y) for x in range(rebox...,nre) for y in range(imbox...,nim)]
    F=map(z->1/f(z),Z); k=isfinite.(F)
    fit=approximate(F[k],Z[k]; max_iter=120,tol=1e-11,stagnation=10)
    Zt=[complex(x,y) for x in range(rebox...,nre-1).+0.007 for y in range(imbox...,nim-1).+0.006]
    res=maximum(abs.(map(fit,Zt).-map(z->1/f(z),Zt)))/maximum(abs.(map(z->1/f(z),Zt)))
    poles(fit), res
end
function arbit(name, f, cands, rebox, imbox)
    println("--- $name ---")
    for (nre,nim) in ((26,22),(30,24),(22,18))
        ps,res=aaa_poles_res(f;rebox,imbox,nre,nim)
        line=@sprintf("   grid %d×%d resid=%.0e:",nre,nim,res)
        for (l,c) in cands
            d,i=findmin(abs.(ps.-c)); line*=@sprintf("  %s Δ=%.3f@(%+.3f%+.3fi)",l,d,real(ps[i]),imag(ps[i]))
        end
        println(line)
    end
end
println("\n== AAA arbiter: closest extrapolated pole to each candidate (fit UHP data) ==")
arbit("mu=2 k=0.5 [Lopez UHP]", ω->ΛL(ω,0.5,2.0),
    [("Ldesc",0.124-0.136im),("VMD",0.16297-0.11048im)],(-0.05,0.45),(0.05,0.55))
arbit("mu=2 k=0.5 [VMD  UHP]", vmd_defl(plasma2,0.5),
    [("Ldesc",0.124-0.136im),("VMD",0.16297-0.11048im)],(-0.05,0.45),(0.05,0.55))
arbit("mu=2 k=0.6 [Lopez UHP]", ω->ΛL(ω,0.6,2.0),
    [("Ldesc",0.004-0.234im),("VMD",0.20540-0.18512im)],(-0.05,0.45),(0.05,0.55))
arbit("mu=2 k=0.6 [VMD  UHP]", vmd_defl(plasma2,0.6),
    [("Ldesc",0.004-0.234im),("VMD",0.20540-0.18512im)],(-0.05,0.45),(0.05,0.55))
# mu=10 descent: Lopez root 0.277-0.988i. VMD root? find it.
g=vmd_defl(plasma10,2.15)
for s in (0.4-0.5im,0.277-0.988im,0.0-1.0im)
    r=muller(g,s-1e-3,s,s+1e-3im); @printf("  μ=10 k=2.15 VMD seed %+.2f%+.2fi -> %+.4f%+.4fi\n",real(s),imag(s),real(r),imag(r))
end
arbit("mu=10 k=2.15 [Lopez UHP]", ω->ΛL(ω,2.15,10.0),
    [("Ldesc",0.277-0.988im),("VMDrise",0.55-0.30im)],(0.0,0.75),(0.1,0.85))
arbit("mu=10 k=2.15 [VMD  UHP]", vmd_defl(plasma10,2.15),
    [("Ldesc",0.277-0.988im),("VMDrise",0.55-0.30im)],(0.0,0.75),(0.1,0.85))
