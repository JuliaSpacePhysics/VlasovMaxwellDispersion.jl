include("lopez.jl")
using VlasovMaxwellDispersion
using VlasovMaxwellDispersion: DispersionFunction
using RationalFunctionApproximation: approximate, poles
using Printf

pair(vdf) = (NormalizedSpecies(1.0, 1.0, vdf), NormalizedSpecies(-1.0, 1.0, vdf))
plasma2 = pair(MaxwellJuttner(mu = 2.0))
plasma10 = pair(MaxwellJuttner(mu = 10.0))

vmd_det(plasma, kz) = DispersionFunction(plasma, Wavenumber(0.0, kz); mode = :L, deflate = false)
vmd_defl(plasma, kz) = DispersionFunction(plasma, Wavenumber(0.0, kz); mode = :L, deflate = true)

println("== VMD roots via Muller on L-mode factor, seeded from both candidates ==")
for (kz, seeds) in ((0.5, [("Lopez-desc",0.124-0.136im),("VMD-riser",0.163-0.110im)]),
                    (0.6, [("Lopez-desc",0.004-0.234im),("VMD-riser",0.20-0.16im)]))
    g = vmd_det(plasma2, kz)
    for (lbl,s) in seeds
        r = muller(g, s-1e-3, s, s+1e-3im)
        @printf("  μ=2 k=%.2f %-11s seed=%+.3f%+.3fim -> VMD root=%+.5f%+.5fim |det|rel=%.1e\n",
            kz, lbl, real(s),imag(s), real(r),imag(r), abs(g(r))/abs(g(s+0.3im)))
    end
end

# --- AAA arbiter ---
# Fit 1/f on a dense UHP grid; poles of the fit = zeros of f (extrapolated to LHP).
function aaa_roots(f; rebox=(-0.3,0.5), imbox=(0.05,0.6), nre=26, nim=22)
    Z = ComplexF64[]
    for xr in range(rebox...,nre), yi in range(imbox...,nim)
        push!(Z, complex(xr,yi))
    end
    F = map(z->1/f(z), Z)
    keep = isfinite.(F)
    fit = approximate(F[keep], Z[keep]; max_iter=120, tol=1e-11, stagnation=10)
    ps = poles(fit)
    # held-out residual: fit on even grid, test on shifted points
    Zt = [complex(xr,yi) for xr in range(rebox...,nre-1).+0.007 for yi in range(imbox...,nim-1).+0.006]
    ft = map(fit, Zt); ref = map(z->1/f(z), Zt)
    res = maximum(abs.(ft.-ref))/maximum(abs.(ref))
    ps, res, fit
end

println("\n== AAA arbiter: fit TRUE UHP data, extrapolate poles to LHP ==")
for (name, kz, cands) in (
        ("mu=2 k=0.5", 0.5, [("Lopez-desc",0.124-0.136im),("VMD-riser",0.163-0.110im)]),
    )
    println("--- $name; candidate LHP roots: ", cands, " ---")
    for (src, f) in (("Lopez ΛL", ω->ΛL(ω,kz,2.0)), ("VMD L-mode", vmd_defl(plasma2,kz)))
        ps, res, _ = aaa_roots(f)
        lhp = filter(p-> -0.6<imag(p)<-0.02 && -0.4<real(p)<0.6, ps)
        @printf("  [%s] held-out resid=%.1e ; LHP poles:\n", src, res)
        for p in sort(lhp, by=imag, rev=true)
            dists = [(l,abs(p-c)) for (l,c) in cands]
            best = argmin(last.(dists))
            @printf("      %+.4f%+.4fim  (nearest %s, Δ=%.3f)\n", real(p),imag(p), dists[best][1], dists[best][2])
        end
    end
end
