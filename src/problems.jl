# SciML-style problem/algorithm/solution layer over the root problem det(𝒟)=0.

abstract type AbstractDispersionProblem end

"""
    LocalDispersionProblem(plasma, k, omega0; closure=HarmonicSum())

Refine the single mode of `det(𝒟(plasma,ω,k))=0` near seed `omega0`.
"""
struct LocalDispersionProblem{P,K,T,C} <: AbstractDispersionProblem
    plasma::P
    k::K
    omega0::T
    closure::C
end
LocalDispersionProblem(plasma, k, omega0; closure=HarmonicSum()) =
    LocalDispersionProblem(plasma, k, complex(float(omega0)), closure)

"""
    GlobalDispersionProblem(plasma, k, region; closure=HarmonicSum())

Find *all* roots and poles of `det(𝒟(plasma,ω,k))` in the complex box
`region=(lowerleft, upperright)`.
"""
struct GlobalDispersionProblem{P,K,R,C} <: AbstractDispersionProblem
    plasma::P
    k::K
    region::R
    closure::C
end
GlobalDispersionProblem(plasma, k, region; closure=HarmonicSum()) =
    GlobalDispersionProblem(plasma, k, region, closure)

"""
    BranchProblem(plasma, ks, omega0; closure=HarmonicSum())

Track one dispersion branch across the wavenumber sequence `ks`, seeded at
`omega0`.`
"""
struct BranchProblem{P,K,T,C} <: AbstractDispersionProblem
    plasma::P
    ks::K
    omega0::T
    closure::C
end
BranchProblem(plasma, ks, omega0; closure=HarmonicSum()) =
    BranchProblem(plasma, ks, complex(float(omega0)), closure)


abstract type DispersionAlgorithm end


"""
    GRPF(; tol=1e-3, params=nothing)

Global Roots-and-Poles Finder (argument principle) for
[`GlobalDispersionProblem`](@ref). `params::GRPFParams` overrides the defaults.
"""
Base.@kwdef struct GRPF <: DispersionAlgorithm
    tol::Float64 = 1.0e-3
    params::Union{Nothing, GRPFParams} = nothing
end

"""
    ArcLength(; atol=1e-10, maxiter=100, fallback=true)

Arc-length branch continuation for [`BranchProblem`](@ref): predictor-seeded
Muller per `k`, with optional local-GRPF fallback on failed/large jumps. See
[`track`](@ref) for the `fallback` NamedTuple knobs.
"""
Base.@kwdef struct ArcLength <: DispersionAlgorithm
    atol::Float64 = 1.0e-10
    maxiter::Int = 100
    fallback = true
end

"""
    DispersionSolution

Result of `solve(::AbstractDispersionProblem, alg)`. `omega` is the root
(`Local`), the vector of roots (`Branch`), or all roots (`Global`); `poles` is
non-`nothing` only for `Global`. `retcode` is `:Success`, `:Failure`, or
`:Partial` (branch with some non-converged `k`). `resid` is the
scale-invariant [`residual`](@ref) `|det 𝒟| / ∏ᵢ‖𝒟ᵢ,:‖` at the root(s),
mirroring the shape of `omega` (`NaN` for non-converged entries). 

Polished roots (`Local`/`Branch`) reach ~machine epsilon; 
`Global` roots are only mesh-accurate so theirs sit near `alg.tol`.
"""
struct DispersionSolution{T, Po, R, Pr, A}
    omega::T
    poles::Po
    resid::R
    retcode::Symbol
    prob::Pr
    alg::A
end

function Base.show(io::IO, sol::DispersionSolution)
    print(io, "DispersionSolution(retcode=:", sol.retcode, ", omega=", sol.omega)
    isnothing(sol.poles) || print(io, ", ", length(sol.poles), " poles")
    print(io, ")")
    return
end
