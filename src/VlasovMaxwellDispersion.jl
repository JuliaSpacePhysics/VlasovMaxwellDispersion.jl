module VlasovMaxwellDispersion

include("../lib/PlasmaBase/src/PlasmaBase.jl")
using .PlasmaBase
using .PlasmaBase: AbstractVDF, AbstractPlasma, Particle, Species, Plasma,
    charge, mass, particle, number_density, distribution, species, magnetic_field, frequency,
    gyrofrequency_ratio, plasma_gyro_ratio
using LinearAlgebra
using StaticArrays
using Bumper: @no_escape, @alloc
using Roots
using RootsAndPoles
using SpecialFunctions
using SpecialFunctions: erfcx, gamma
using QuadGK
using NonNegLeastSquares: nonneg_lsq
using CommonSolve
import CommonSolve: solve

function contribution end

include("derivatives.jl")
include("holomorphic_ad.jl")   # erfc/erf/gamma/… differentiable at complex argument
include("Bessel.jl")
include("types.jl")
include("problems.jl")
include("integrals.jl")
include("hilbert_pwpoly.jl")   # parallel H∥ piecewise-poly primitive
include("perp_analytic.jl")    # perpendicular P⊥ Bessel-moment primitive
include("projection.jl")
include("distributions/distributions.jl")
include("builders.jl")         # particle-identity + physical-unit Species adapters
include("susceptibility.jl")
include("solve.jl")
include("track.jl")

export Regime, NonRelativistic, Relativistic
export Wavenumber, para, perp
# re-exported physical vocabulary from PlasmaBase
export AbstractVDF, Particle, Electron, Proton, Ion, Species, Plasma
export NormalizedSpecies
export Maxwellian, MaxwellJuttner, ColdVDF, GridVDF, SeparableVDF, ReducedVDF, CoupledVDF, GaussianRing
export Separable, ⊗, Gaussian, GyroRing
export plasma_dispersion_function, Z, hilbert
export IntegralClosure, HarmonicSum, Newberger
export contribution, dielectric, dispersion_tensor, 𝒟, electrostatic_det, solve, residual, dispersion_residual
export LocalDispersionProblem, GlobalDispersionProblem, BranchProblem
export DispersionAlgorithm, Muller, Secant, GRPF, ArcLength, DispersionSolution
export GridFitMethod, NonnegBSpline, BicubicHermite, fit_grid

end
