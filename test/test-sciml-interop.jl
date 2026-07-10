# det 𝒟 is holomorphic but complex-valued, so most NonlinearSolve solvers break: 
# the bracketing methods need a real sign change except `BracketingNonlinearSolve.Muller`
# and common SimpleNonlinearSolve method Broyden/DFSane compare residual magnitudes via `isless(::Complex,…)`
# except SimpleHalley (FiniteDiff AD).

using Test
using VlasovMaxwellDispersion

import SimpleNonlinearSolve as SNS
import BracketingNonlinearSolve as BNS

p = NormalizedSpecies(Omega=1.0, Pi2=1.0, vdf=Maxwellian(1.0))
k = Wavenumber(kperp=0.01, kz=0.5)
prob = DispersionProblem(p, 0.6, k)
f = prob.f
ref = solve(prob).omega
@test abs(f(ref)) < 1e-8

h = 1e-3
# IntervalNonlinearProblem/Muller; seed the imaginary middle so
# it tracks the complex root (default middle (l+r)/2 is real).
ip = SNS.IntervalNonlinearProblem((ω, _) -> f(ω), (prob.omega0 - h, prob.omega0 + h))
sm = solve(ip, BNS.Muller(prob.omega0 + h * im))
@test sm.retcode == SNS.ReturnCode.Success
@test sm.u ≈ ref rtol=1e-6

# Initial-guess path
np = SNS.NonlinearProblem((ω, _) -> f(ω), prob.omega0)
sh = solve(np, SNS.SimpleHalley())
@test sh.u ≈ ref rtol=1e-6

# SciMLBaseExt with SciML alg
sd = solve(prob, SNS.SimpleHalley())
@test sd.retcode == :Success
@test sd.omega ≈ ref rtol=1e-6
