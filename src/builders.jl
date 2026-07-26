# Bridge: a PHYSICAL Plasma → this solver's dimensionless rep

import .PlasmaBase: scales

scales(p::NormalizedPlasma) = p.scales

"""
    NormalizedPlasma(plasma::Plasma)

Normalize a physical plasma into this solver's dimensionless form: `Ω̃`, `Π̃²` per species
from `B0` and the plasma's `ref`, and every VDF spec resolved against the resulting
[`scales`](@ref).
"""
function NormalizedPlasma(plasma::Plasma)
    sps = Tuple(species(plasma))
    sc = scales(plasma)
    return NormalizedPlasma(ntuple(i -> NormalizedSpecies(sps[i], scales(sc, i)), length(sps)), sc)
end

prepare(p::Plasma, args...; kw...) = prepare(NormalizedPlasma(p), args...; kw...)
