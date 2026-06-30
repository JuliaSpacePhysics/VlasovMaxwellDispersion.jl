# Perf experiment for nonrelativistic CoupledVDF harmonic path
# Run: julia --project=. benchmark/coupled_exp.jl
#

# , then exploit
# that `_In_block` is linear in M=c·Δm → scale the 3 momentum moments (3 muls)
# instead of premultiplying the 5-vector (5 muls). See `_coupled_block`.
# Tried: accumulating the 6 distinct tensor scalars instead of summing SMatrices —
# the chained `a=b=…=0` init boxed inside the quadgk closure (1.7M allocs); reverted.
# Tried: factoring real Bessel weights before `im*…` — slightly slower; reverted.
using VlasovMaxwellDispersion
using VlasovMaxwellDispersion: contribution
using Chairmarks, Statistics

g0(w, u) = exp(-(u^2 + w^2 + 0.6u * w))
s = NormalizedSpecies(-1.0, 1.0, CoupledVDF(g0; para = (-8.0, 8.0), perp = 6.0))
ω = 1.2 + 0.05im
for kp in (0.3, 1.0, 2.0)
    k = Wavenumber(kp, 0.4)
    contribution(s, ω, k)   # warm
    b = @be contribution($s, $ω, $k) seconds = 2
    println("kperp=$kp  min=$(round(minimum(b).time * 1e3, digits = 3))ms  ",
        "median=$(round(median(b).time * 1e3, digits = 3))ms  allocs=$(minimum(b).allocs)")
end
