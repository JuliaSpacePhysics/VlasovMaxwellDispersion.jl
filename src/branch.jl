"""
Dispersion branch for parameter sweep. 
For a fixed wavevector, `omega`, `k`, and `resid` are scalars. 
Otherwise they are arrays shaped like the parameter grid.
Missing branch points are `NaN` in `omega` and `resid`.
"""
struct DispersionBranch{W,K,R}
    omega::W
    k::K
    resid::R
end

Base.length(b::DispersionBranch) = length(b.omega)
Base.iterate(b::DispersionBranch, args...) = iterate(b.omega, args...)
Base.getindex(b::DispersionBranch, args...) = getindex(b.omega, args...)
Base.eltype(::Type{<:DispersionBranch{W}}) where {W} = eltype(W)
Base.lastindex(b::DispersionBranch) = lastindex(b.omega)

isgrowing(b, margin = 1e-3) = maximum((imag(ω) for ω in b if isfinite(ω)); init=(-Inf)) > margin
