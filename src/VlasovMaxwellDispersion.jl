module VlasovMaxwellDispersion

using PlasmaPopulations
using PlasmaPopulations: C_SI, si_frequency, kappa_widths, ParameterGeometry, gyrofrequency, plasma_frequency, alfven_speed
using LinearAlgebra
using StaticArrays
using Bumper: @no_escape, @alloc
using Gamma: gamma
import Gamma
using SpecialFunctions
using SpecialFunctions: erfcx
import SpecialFunctions as SF
using HypergeometricFunctions: _₂F₁
using QuadGK
using NonNegLeastSquares: nonneg_lsq
using CommonSolve
import CommonSolve: solve, init, solve!, step!

include("api.jl")
include("utils.jl")
include("derivatives.jl")
include("erased.jl")
include("Bessel.jl")
include("types.jl")
include("branch.jl")
include("reduction.jl")
include("problems.jl")
include("integrals.jl")
include("perp_analytic.jl")    # perpendicular P⊥ Bessel-moment primitive
include("projection.jl")
include("table_plan.jl")
include("distributions/distributions.jl")
include("builders.jl")         # particle-identity + physical-unit Population adapters
include("assemble.jl")
include("susceptibility.jl")
include("dispersion_function.jl")
include("CoupledPlan.jl")
include("solve.jl")

export Particle, Electron, Positron, Proton, Population, Plasma,
    BiMaxwellian, GaussianRing, BiKappa, ProductBiKappa, Cold, MaxwellJuttner, GridVDF,
    Wavevector, AngleSweep, CartesianSweep
export NonRelativistic, Relativistic
export NormalizedPopulation, NormalizedPlasma, scales
export SplineVDF, SeparableVDF, ReducedVDF, CoupledVDF, LowRankVDF
export ⊗, Gaussian, GyroRing, Kappa
export HarmonicSum, Newberger, prepare
export FixedNodeEval, AdaptiveEval
export contribution, dispersion_tensor, electrostatic_det, residual
export DispersionProblem, Seed, Region, ReturnCode, solve
export Muller, GRPF, Continuation, AAA
export isgrowing, dispersion_diagram!, dispersion_diagram
export NonnegBSpline, BicubicHermite

VERSION >= v"1.11.0-DEV.469" && eval(Meta.parse("""public Regime, SolverVDF, IntegralClosure, ChiBackend, Separable, TensorReduction,
    𝒟, plan_contribution, fit_grid, DispersionFunction,
    DispersionSolution, SurveySolution, DispersionBranch, SolveStats, successful_retcode, init, solve!, step!"""))

function dispersion_diagram! end
function dispersion_diagram end

include("precompile.jl")

end
