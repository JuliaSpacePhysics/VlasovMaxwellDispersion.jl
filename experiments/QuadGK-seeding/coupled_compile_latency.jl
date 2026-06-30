# Cold compile-latency probe for the nonrelativistic CoupledVDF harmonic path.
# Measures `@timed(...).compile_time` of a first coupled `contribution` in a
# FRESH process (compilation isolated from runtime). Measures whatever is checked
# out: ~3–4 s/case on main; `git apply coupled_seeding.patch` then rerun for the
# seeded ~15–17 s/case. See coupled_compile_latency.md for the analysis.

const PROJ = dirname(Base.active_project())

# One representative coupled case per distinct seed path. Each runs in its own
# process — a shared process would amortize the very compilation we measure.
const CASES = [
    ("inseparable kp=0.3",
        "g0(u,v)=exp(-(u^2+v^2+0.6u*v)); s=NormalizedSpecies(-1.0,1.0,CoupledVDF(g0;para=(-8.0,8.0),perp=6.0)); contribution(s,1.2+0.05im,Wavenumber(0.3,0.4))"),
    ("perp-cut ring kp=0.6",
        "f0=Maxwellian(vth_para=0.1,vth_perp=0.05,vr=0.6); s=CoupledVDF(f0;para=(-0.8,0.8),perp=(0.3,1.05)); contribution(s,1.3+0.02im,Wavenumber(0.6,0.4))"),
]

probe(call) = """
using VlasovMaxwellDispersion
using VlasovMaxwellDispersion: CoupledVDF, NormalizedSpecies, Maxwellian, Wavenumber, contribution
r = @timed ($call)
print(r.compile_time)
"""

println("cold compile_time (fresh process per case):")
for (name, call) in CASES
    out = read(`$(Base.julia_cmd()) --project=$PROJ -e $(probe(call))`, String)
    println("  ", rpad(name, 22), round(parse(Float64, out), digits = 3), " s")
end
