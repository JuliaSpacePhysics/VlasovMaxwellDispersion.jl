using RootsAndPoles

"""
    GRPF(; tol=1e-3, meshtol=nothing, params=nothing)

Global Roots-and-Poles Finder (argument principle) for [`GlobalDispersionProblem`](@ref).
`tol` is the refinement accuracy; `meshtol` is the *initial* mesh spacing GRPF refines from
(`nothing` ⇒ a coarse fraction of the box). `params::GRPFParams` overrides the defaults.

Works on the deflated `det(ω̃²𝒟)`, so genuine roots within `2tol` of `ω=0` are
dropped together with the deflation's origin artifact.
"""
Base.@kwdef struct GRPF{P}
    tol::Float64 = 1.0e-3
    meshtol::Union{Nothing, Float64} = nothing
    params::P = nothing
end

struct GRPFCache{P, A, K, F, R}
    prob::P
    alg::A
    k::K
    f::F
    nevals::Base.RefValue{Int}
    refine::R
end

function CommonSolve.init(prob::GlobalDispersionProblem, alg::GRPF; refine = Muller())
    iszero(n_swparams(prob)) ||
        throw(ArgumentError("GRPF solves only the fixed-k (m=0) problem; use DualAAA for parameter sweeps"))
    k_fixed = wavefun(prob)()
    nev = Ref(0)
    f = ω -> (nev[] += 1; det(wave_dispersion_tensor(prob.plasma, ω, k_fixed; closure = prob.closure)))
    return GRPFCache(prob, alg, k_fixed, f, nev, refine)
end

function CommonSolve.solve!(cache::GRPFCache)
    (; prob, alg, f, refine) = cache
    roots, poles = _grpf_roots(f, prob.region; alg.tol, alg.meshtol, alg.params)
    # Drop the deflated det's spurious ω=0 zero (perp k): it sits at |ω| ≲ tol
    # (mesh accuracy); 2·tol gives a one-cell margin.
    _in_box(prob.region) && filter!(ω -> abs(ω) > 2alg.tol, roots)
    retcode = isempty(roots) ? :Failure : :Success
    return _fixedk_survey(prob, alg, cache.k, roots, poles, cache.nevals[], retcode; refine)
end

# ComplexF64 throughout: RootsAndPoles' Delaunay geometry (IndexablePoint2D) is hard Float64
function _grpf_roots(f, region; tol = 1.0e-3, meshtol = nothing, params = nothing)
    lowerleft, upperright = ComplexF64(region[1]), ComplexF64(region[2])
    diag = hypot(real(upperright) - real(lowerleft), imag(upperright) - imag(lowerleft))
    base = @something meshtol max(diag / 24, tol)
    origcoords = rectangulardomain(lowerleft, upperright, base)
    p = @something params GRPFParams(5000, tol, false)
    return grpf(f, origcoords, p)
end
