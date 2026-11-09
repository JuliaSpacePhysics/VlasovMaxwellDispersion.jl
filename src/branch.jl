"""
Dispersion branch for parameter sweep. 
For a fixed wavevector, `omega`, `k`, and `resid` are scalars. 
Otherwise they are arrays shaped like the parameter grid.
Missing branch points are `NaN` in `omega` and `resid`. 
`omega` and `k` are in the units of the problem's plasma; `omega_unit` is `Ω_ref` in those units.
"""
struct DispersionBranch{W,K,R,U}
    omega::W
    k::K
    resid::R
    omega_unit::U
end

Base.length(b::DispersionBranch) = length(b.omega)
Base.iterate(b::DispersionBranch, args...) = iterate(b.omega, args...)
Base.getindex(b::DispersionBranch, args...) = getindex(b.omega, args...)
Base.eltype(::Type{<:DispersionBranch{W}}) where {W} = eltype(W)
Base.lastindex(b::DispersionBranch) = lastindex(b.omega)

isgrowing(b::DispersionBranch, margin = 1e-3) =
    maximum((imag(ω) for ω in b if isfinite(ω)); init=(-Inf)) > margin * b.omega_unit
