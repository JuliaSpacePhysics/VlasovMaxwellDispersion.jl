# The observable: field-aligned current needed to destabilise the ion-cyclotron band as a
# function of how transversely heated the O⁺ already is. u_c → 0 is the collapse.

include(joinpath(@__DIR__, "common.jl"))

const VRS = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5)
const DELTAS = (0.1, 0.3, 1.0)

rows = NTuple{9,Any}[]
for δ in DELTAS, vr in VRS
    uc, g, kr, r, wr = threshold(u -> plasma(; delta = δ, vr, u))
    j = isnan(uc) ? NaN : current_density(uc) * 1e6
    push!(rows, (δ, vr, ring_energy(vr), uc, uc * sqrt(MI_ME), j, kr, r, wr))
    @printf("δ=%.1f v_r/w=%.1f E⊥=%5.2f eV  u_c=%.4f v_e (%5.1f w_i)  j_c=%7.3f µA/m²  ω_r=%.3f\n",
        δ, vr, ring_energy(vr), uc, uc * sqrt(MI_ME), j, wr)
end
writetsv(joinpath(@__DIR__, "out_joint.tsv"),
    "delta\tvr\tE_perp_eV\tu_c\tu_c_wi\tj_uA\tkrho\tratio\tomega_r", rows)
