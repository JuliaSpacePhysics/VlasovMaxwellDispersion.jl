"""
    link(vals; gate)

Link per-grid-point complex value sets `vals[I]` into continuous sheets over
the m-D parameter grid.
Values match extrapolated tips of sheets at neighboring points greedily within
`gate`; unmatched values start new sheets.
"""
function link(vals::AbstractArray{<:Any, M}; gate) where {M}
    V = ComplexF64
    Tip = Tuple{Int, V, V}                       # (sheet id, value, velocity)
    sheets = Vector{Tuple{CartesianIndex{M}, V}}[]
    tips = Dict{CartesianIndex{M}, Vector{Tip}}()
    for c in CartesianIndices(size(vals))
        vs = V.(vals[c])
        # tip of each sheet at the immediate predecessor along every axis
        ref = Dict{Int, Tuple{V, V}}()           # sid ⇒ (tip, velocity)
        for j in 1:M
            nb = c - CartesianIndex(ntuple(i -> Int(i == j), M))
            nb[j] ≥ 1 || continue
            for (sid, v, dv) in get(tips, nb, Tip[])
                haskey(ref, sid) || (ref[sid] = (v, dv))
            end
        end
        # greedy: smallest |value − extrapolated tip| first, one value per sheet
        cand = sort!(
            [
                (abs(vs[i] - (v + dv)), i, sid, v)
                    for (sid, (v, dv)) in ref for i in eachindex(vs)
            ]; by = first
        )
        taken = falses(length(vs))
        assigned = Tip[]
        for (d, i, sid, v) in cand
            (d ≤ gate && !taken[i] && all(a -> a[1] ≠ sid, assigned)) || continue
            taken[i] = true
            push!(sheets[sid], (c, vs[i]))
            push!(assigned, (sid, vs[i], vs[i] - v))
        end
        for i in eachindex(vs)                   # leftovers start new sheets
            taken[i] && continue
            push!(sheets, [(c, vs[i])])
            push!(assigned, (length(sheets), vs[i], zero(V)))
        end
        tips[c] = assigned
    end
    return map(sheets) do sh
        out = fill(complex(NaN, NaN), size(vals))
        foreach(((c, v),) -> out[c] = v, sh)
        out
    end
end
