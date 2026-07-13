"""
    link(vals; gate)

Link per-grid-point complex value sets `vals[I]` into continuous sheets over the
m-D parameter grid.

One sheet is traced at a time, seeded at the highest-growth-rate 
unclaimed value and grown outward: each frontier cell takes the
nearest unclaimed value within `gate` of a quadratic extrapolation.
"""
function link(vals::AbstractArray{<:Any,M}; gate) where {M}
    V = ComplexF64
    grid = CartesianIndices(size(vals))
    pool = map(c -> V.(vals[c]), grid)               # candidate values per cell
    live = map(p -> trues(length(p)), pool)          # still unclaimed
    steps = [CartesianIndex(ntuple(i -> ifelse(i == j, s, 0), M)) for j in 1:M for s in (-1, 1)]
    sheets = Vector{Pair{CartesianIndex{M},V}}[]
    while true
        sc = zero(CartesianIndex{M})                 # seed at the highest-growth live value
        si = 0
        best = -Inf
        for c in grid, i in eachindex(pool[c])
            live[c][i] && imag(pool[c][i]) > best && (best=imag(pool[c][i]); sc=c; si=i)
        end
        si == 0 && break
        live[sc][si] = false
        claim = Dict(sc => pool[sc][si])
        frontier = [sc + s for s in steps if (sc + s) in grid]
        while !isempty(frontier)
            c = popfirst!(frontier)
            (c in grid && !haskey(claim, c)) || continue
            w = _quadratic_extrapolate(claim, c, steps)
            i = isnothing(w) ? 0 : _nearest_live(pool[c], live[c], w, gate)
            i == 0 && continue
            live[c][i] = false
            claim[c] = pool[c][i]
            for s in steps
                nb = c + s
                (nb in grid && !haskey(claim, nb)) && push!(frontier, nb)
            end
        end
        push!(sheets, collect(claim))
    end
    return map(sheets) do sh
        out = fill(complex(NaN, NaN), size(vals))
        foreach(((c, v),) -> out[c] = v, sh)
        out
    end
end

# Quadratic extrapolation lets a smoothly bending branch cross a flat one without swapping
# when the bend is mild. A kink where the second difference `|a−2b+d|` exceeds `_CURV·|a−b|` 
# (e.g. an instability bifurcation) would overshoot onto a neighbouring branch -> linear
const _CURV = 0.4
function _quadratic_extrapolate(claim, c, steps)
    acc = zero(valtype(claim))
    n = 0
    for s in steps
        a = get(claim, c + s, nothing)
        isnothing(a) && continue
        b = get(claim, c + 2s, nothing)
        d = get(claim, c + 3s, nothing)
        acc += if isnothing(b)
            a
        elseif isnothing(d) || abs(a - 2b + d) > _CURV * abs(a - b)
            2a - b
        else
            3a - 3b + d
        end
        n += 1
    end
    return n == 0 ? nothing : acc / n
end

# Index of the nearest unclaimed value to `w` within `gate`; 0 if none.
function _nearest_live(vs, live, w, gate)
    best = 0
    bd = gate
    for i in eachindex(vs)
        live[i] || continue
        d = abs(vs[i] - w)
        d ≤ bd && (bd=d; best=i)
    end
    return best
end
