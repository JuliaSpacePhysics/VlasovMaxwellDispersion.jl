# Currentless growth: how sharp must the O⁺ ring be before it destabilises the
# ion-cyclotron band on its own? Electrons carry no drift anywhere in this file.

include(joinpath(@__DIR__, "common.jl"))

rows = map(0.0:0.25:3.5) do vr
    g, kr, r, wr = gamma_max(plasma(; vr))
    ok = confirm(plasma(; vr), g, kr, r, wr)
    (vr, ring_energy(vr), g, kr, r, wr, ok)
end

println("currentless ring-driven growth (u = 0), full ring population δ = 1")
println("  v_r/w   E⊥[eV]     γ/Ω     k⊥ρ   k∥/k⊥    ω_r/Ω   confirmed")
for (vr, E, g, kr, r, wr, ok) in rows
    @printf("%7.2f  %7.2f  %+8.5f  %5.2f  %6.4f  %7.3f   %s\n", vr, E, g, kr, r, wr, ok)
end
writetsv(joinpath(@__DIR__, "out_ring.tsv"),
    "vr\tE_perp_eV\tgamma\tkrho\tratio\tomega_r\tconfirmed", rows)

vc = threshold(vr -> plasma(; vr); lo = 0.0, hi = 3.5, iters = 12)
@printf("\nthreshold ring speed: v_r/w = %.3f  →  E⊥ = %.2f eV  (γ = 1e-3 Ω, ω_r = %.3f Ω)\n",
    vc[1], ring_energy(vc[1]), vc[5])

# How much of the population has to be ring-like?
println("\nring fraction δ needed at fixed v_r (currentless):")
for vr in (2.0, 2.5, 3.0)
    dc = threshold(δ -> plasma(; delta = δ, vr); lo = 0.02, hi = 1.0, iters = 10)
    @printf("  v_r/w = %.1f (E⊥ = %.1f eV):  δ_c = %.3f\n", vr, ring_energy(vr), dc[1])
end
