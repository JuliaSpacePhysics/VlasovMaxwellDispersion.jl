module VlasovMaxwellDispersion

include("../lib/PlasmaBase/src/PlasmaBase.jl")
using .PlasmaBase
using .PlasmaBase: AbstractVDF, AbstractPlasma, Particle, Species, Plasma,
    charge, mass, particle, number_density, distribution, species, magnetic_field, frequency,
    gyrofrequency_ratio, plasma_gyro_ratio
using LinearAlgebra
using StaticArrays
using Bumper: @no_escape, @alloc
using SpecialFunctions
using SpecialFunctions: erfcx, gamma
using HypergeometricFunctions: _₂F₁
using QuadGK
using NonNegLeastSquares: nonneg_lsq
using CommonSolve
import CommonSolve: solve, init, solve!, step!

function contribution end

include("derivatives.jl")
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

export Regime, NonRelativistic, Relativistic
export Wavenumber, para, perp
# re-exported physical vocabulary from PlasmaBase
export AbstractVDF, Particle, Electron, Proton, Species, Plasma
export NormalizedSpecies
export Maxwellian, MaxwellJuttner, ColdVDF, GridVDF, SeparableVDF, ReducedVDF, CoupledVDF, GaussianRing
export BiKappa, ProductBiKappa, Kappa
export Separable, ⊗, Gaussian, GyroRing
export plasma_dispersion_function, Z
export IntegralClosure, HarmonicSum, Newberger, prepare
export contribution, dielectric, dispersion_tensor, 𝒟, electrostatic_det, solve, init, solve!, step!, residual
export DispersionProblem, GlobalDispersionProblem
export Muller, GRPF, ArcLength, AAA, JumpFallback, DispersionSolution
export AngleSweep, CartesianSweep
export SurveySolution, DispersionBranch, SolveStats, dispersion_diagram
export GridFitMethod, NonnegBSpline, BicubicHermite, fit_grid

function dispersion_diagram end

end
