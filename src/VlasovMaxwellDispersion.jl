module VlasovMaxwellDispersion

using LinearAlgebra
using StaticArrays
using Roots
using RootsAndPoles
using SpecialFunctions
using SpecialFunctions: erfcx, gamma
using QuadGK, HCubature
using NonNegLeastSquares
using CommonSolve
import CommonSolve: solve
using ForwardDiff

function contribution end

include("derivatives.jl")
include("Bessel.jl")
include("types.jl")
include("problems.jl")
include("integrals.jl")
include("hilbert_pwpoly.jl")   # parallel H∥ piecewise-poly primitive
include("perp_analytic.jl")    # perpendicular P⊥ Bessel-moment primitive
include("projection.jl")
include("distributions.jl")
include("susceptibility.jl")
include("solve.jl")
include("track.jl")

export Regime, NonRelativistic, Relativistic, Continuation, Separability
export Wavenumber, para, perp, Species, Plasma
export Maxwellian, MaxwellJuttner, ColdVDF, GridVDF, SeparableVDF, CoupledVDF
export plasma_dispersion_function, Z, Gamma_n, hilbert
export IntegralClosure, HarmonicSum, Newberger, besselj_complex
export contribution, dielectric, dispersion_tensor, 𝒟, electrostatic_det, solve, residual
export LocalDispersionProblem, GlobalDispersionProblem, BranchProblem
export DispersionAlgorithm, Muller, Secant, GRPF, ArcLength, DispersionSolution
export GridFitMethod, NonnegBSpline, BicubicHermite, fit_grid

end
