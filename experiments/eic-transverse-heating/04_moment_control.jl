# Control: is the ring's effect visible in the moments an observer reports? Each row pairs
# a ring with a bi-Maxwellian of identical n, T∥ and T⊥ = (1+v_r²)T∥ — same first two
# velocity moments, different shape.

include(joinpath(@__DIR__, "common.jl"))

rows = map((0.0, 0.5, 1.0, 1.5, 2.0, 2.5)) do vr
    A = 1 + vr^2
    gr = gamma_max(plasma(; vr))
    gm = gamma_max(moment_matched(; vr))
    ur = threshold(u -> plasma(; vr, u))[1]
    um = threshold(u -> moment_matched(; vr, u))[1]
    @printf("v_r/w=%.1f  T⊥/T∥=%5.2f | ring: γ=%+.4f u_c=%.4f | bi-Max: γ=%+.4f u_c=%.4f\n",
        vr, A, gr[1], ur, gm[1], um)
    (vr, A, gr[1], gr[4], ur, gm[1], gm[4], um)
end
writetsv(joinpath(@__DIR__, "out_control.tsv"),
    "vr\tA\tgamma_ring\tomega_ring\tu_c_ring\tgamma_bimax\tomega_bimax\tu_c_bimax", rows)
