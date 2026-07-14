using RootsAndPoles

"""
    GRPF(; tol=1e-3, meshtol=nothing, params=nothing)

Global Roots-and-Poles Finder (argument principle) solver. `tol` is the
refinement accuracy; `meshtol` the *initial* mesh spacing (`nothing` ⇒ a coarse
fraction of the box); `params::GRPFParams` overrides the defaults. Searches the
deflated pole-free `det(ω̃²𝒟)` and drops roots within `2tol` of `ω=0` together
with the deflation's origin artifact.
"""
Base.@kwdef struct GRPF{P}
    tol::Float64 = 1.0e-3
    meshtol::Union{Nothing, Float64} = nothing
    params::P = nothing
end

function discover(alg::GRPF, f0, region; keep = Returns(true))
    n = Ref{Int}(0)
    f = ω -> (n[] += 1; f0(ω))
    roots, _ = _grpf_roots(f, region; alg.tol, alg.meshtol, alg.params)
    return roots, n[]
end

# GRPF locates the origin artifact only to mesh accuracy (|ω| ≲ tol even before
# polish); 2tol gives a one-cell margin — genuine roots inside it are dropped too.
_origin_gate(alg::GRPF, diag) = 2 * alg.tol

# ComplexF64 throughout: RootsAndPoles' Delaunay geometry (IndexablePoint2D) is hard Float64
function _grpf_roots(f, region; tol = 1.0e-3, meshtol = nothing, params = nothing)
    lowerleft, upperright = ComplexF64(region[1]), ComplexF64(region[2])
    diag = hypot(real(upperright) - real(lowerleft), imag(upperright) - imag(lowerleft))
    base = @something meshtol max(diag / 24, tol)
    origcoords = rectangulardomain(lowerleft, upperright, base)
    p = @something params GRPFParams(5000, tol, false)
    return grpf(f, origcoords, p)
end
