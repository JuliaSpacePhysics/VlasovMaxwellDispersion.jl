# Gate: reproduce the Kindel & Kennel (1971) current-driven EIC threshold for an O⁺
# topside plasma, and check it is insensitive to the density (kλ_De ≪ 1 here).
#
# Expected from K&K: the unstable root sits just above the first ion gyroharmonic
# (ω_r ≈ 1.2 Ω_i) at k⊥ρ_i ≈ 1-2, with a threshold drift of order 10-20 ion thermal
# speeds at Te/Ti = 1, falling as Te/Ti rises.

include(joinpath(@__DIR__, "common.jl"))

rows = map((0.5, 1.0, 2.0, 4.0, 8.0)) do tau
    uc, g, kr, r, wr = threshold_u(u -> plasma(; u, tau))
    ok = confirm(plasma(; u = uc, tau), g, kr, r, wr)
    (tau, uc, uc * sqrt(tau * MI_ME), current_density(uc; tau) * 1e6, kr, r, wr, g, ok)
end

println("current-driven EIC threshold (γ = 1e-3 Ω_O⁺), pure O⁺, δ=0")
println(" Te/Ti   u_c/v_e   u_c/w_i   j [µA/m²]   k⊥ρ   k∥/k⊥    ω_r     γ      ok")
for (tau, uc, ui, j, kr, r, wr, g, ok) in rows
    @printf("%5.1f  %8.4f  %8.2f  %9.2f  %5.2f  %6.3f  %6.3f  %.4f  %s\n",
        tau, uc, ui, j, kr, r, wr, g, ok)
end
writetsv(joinpath(@__DIR__, "out_gate.tsv"),
    "tau\tu_c\tu_c_wi\tj_uA\tkrho\tratio\tomega_r\tgamma\tconfirmed", rows)

# density insensitivity: same threshold at 10× and 1/10× the reference density
println("\ndensity check at Te/Ti = 1 (reference Π² = $(round(PI2, digits=1)))")
scale_density(s, f) = map(x -> NormalizedSpecies(x.Omega, f * x.Pi2, x.vdf), s)
for f in (0.1, 1.0, 10.0)
    uc = threshold_u(u -> scale_density(plasma(; u), f))[1]
    @printf("  n/n₀ = %4.1f:  u_c/v_e = %.4f\n", f, uc)
end
