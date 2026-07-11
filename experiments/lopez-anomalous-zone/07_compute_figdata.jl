# Figure data for report.typ: CSVs consumed by 08_make_figures.jl
include("lopez.jl")
include("06_corrected_continuation.jl")  # ΛL_corr, dzbar
using VlasovMaxwellDispersion
using VlasovMaxwellDispersion: DispersionFunction
using Printf

outcsv(name, header, rows) = open(joinpath(@__DIR__, name), "w") do io
    println(io, header)
    for r in rows
        println(io, join(r, ","))
    end
end

# --- traces: López-formula root vs VMD root, μ=2, k=0.1→0.85 ---
pair(vdf) = (NormalizedSpecies(1.0, 1.0, vdf), NormalizedSpecies(-1.0, 1.0, vdf))
plasma2 = pair(MaxwellJuttner(mu = 2.0))
vmd_det(kz) = DispersionFunction(plasma2, Wavenumber(0.0, kz); mode = :L, deflate = false)

ks = 0.1:0.05:0.85
function trace(fset, seed)
    out = ComplexF64[]
    s = seed
    for k in ks
        f = fset(k)
        s = muller(f, s - 1.0e-3, s, s + 1.0e-3im)
        push!(out, s)
    end
    return out
end
tl = trace(k -> (ω -> ΛL(ω, k, 2.0)), complex(0.039185, -2.5e-5))
tv = trace(vmd_det, complex(0.039185, -2.5e-5))
outcsv(
    "figdata_traces.csv", "k,lopez_wr,lopez_gm,vmd_wr,vmd_gm",
    [(k, real(a), imag(a), real(b), imag(b)) for (k, a, b) in zip(ks, tl, tv)]
)
println("traces done")

# --- holomorphy maps at k=0.5, μ=2: log10(|∂z̄ f|/|∂z f|) ---
res, ims = range(-0.02, 0.42, 40), range(-0.32, 0.32, 33)
for (name, f) in (("lopez", ω -> ΛL(ω, 0.5, 2.0)), ("corr", ω -> ΛL_corr(ω, 0.5, 2.0)))
    rows = Vector{NTuple{3, Float64}}()
    for x in res, y in ims
        abs(y) < 2.0e-3 && (y = y < 0 ? -2.0e-3 : 2.0e-3)  # dzbar stencil must not straddle the axis
        zb, zz = dzbar(f, complex(x, y); h = 5.0e-4)
        push!(rows, (x, y, log10(max(abs(zb) / abs(zz), 1.0e-9))))
    end
    outcsv("figdata_holo_$name.csv", "re,im,logdefect", rows)
    println("holo $name done")
end
