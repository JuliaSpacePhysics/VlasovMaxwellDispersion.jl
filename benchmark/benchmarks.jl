using BenchmarkTools
using VlasovMaxwellDispersion
using VlasovMaxwellDispersion: contribution
import SimpleNonlinearSolve as SNS
import BracketingNonlinearSolve as BNS
const SUITE = BenchmarkGroup()

# Maxwellian baseline (closed Z/Γ_n harmonic sum) — the fast reference path.
let
    s = Species(-1.0, 1.0, Maxwellian(vth_par=0.9, vth_perp=1.2))
    g = SUITE["maxwellian"] = BenchmarkGroup()
    for kp in (0.1, 1.0)
        k = Wavenumber(kp, 0.4)
        g["kperp=$kp"] = @benchmarkable contribution($s, 1.3 - 0.05im, $k)
    end
end

# CoupledVDF: the two closures (derivation.md §3) on an inseparable f₀.
let
    g0(u, v) = exp(-(u^2 + v^2 + 0.6u * v))
    dpar(u, v) = -(2u + 0.6v) * g0(u, v)
    dperp(u, v) = -(2v + 0.6u) * g0(u, v)
    kw = (parlower=-8.0, parupper=8.0, perpupper=6.0, dpar=dpar, dperp=dperp)
    s = Species(-1.0, 1.0, CoupledVDF(g0; kw...))
    ω = 1.2 + 0.05im
    g = SUITE["coupled_nonrel"] = BenchmarkGroup()
    for kp in (0.3, 1.0)
        k = Wavenumber(kp, 0.4)
        g["B_truncated_kperp=$kp"] = @benchmarkable contribution($s, $ω, $k)
        g["A_newberger_kperp=$kp"] = @benchmarkable contribution($s, $ω, $k; closure=Newberger())
    end
end

let regime=Relativistic()
    μ = 40.0
    γ(u, w) = sqrt(1 + u^2 + w^2)
    f0(u, w) = exp(-μ * γ(u, w))
    dpar(u, w) = -μ * f0(u, w) * u / γ(u, w)
    dperp(u, w) = -μ * f0(u, w) * w / γ(u, w)
    L = sqrt((1 + 25 / μ)^2 - 1)
    kw = (parlower=(-L), parupper=L, perpupper=L, dpar=dpar, dperp=dperp)
    s = Species(1.0, 1.0, CoupledVDF(f0; kw...); regime)
    k = Wavenumber(0.7, 0.4)
    ω = 0.3 + 0.05im
    g = SUITE["Relativistic"] = BenchmarkGroup()
    g["coupled/A_newberger"] = @benchmarkable contribution($s, $ω, $k; closure=Newberger())
    g["coupled/B_truncated"] = @benchmarkable contribution($s, $ω, $k)

    ppar = collect(range(-L, L, length=81))
    pperp = collect(range(0.0, L, length=61))
    F = [f0(u, w) for u in ppar, w in pperp]
    s = Species(1.0, 1.0, GridVDF(ppar, pperp, F; tol=1e-4); regime)
    g["gridvdf"] = @benchmarkable contribution($s, $ω, $k)
end

let
    s = Species(-1.0, 1.0, Maxwellian(vth_par=0.9, vth_perp=1.2))
    k = Wavenumber(0.01, 0.5)
    prob = LocalDispersionProblem(Plasma(s), k, 0.6)
    f = residual(prob)
    h = 1e-3
    g = SUITE["local_solve"] = BenchmarkGroup()
    g["muller_native"] = @benchmarkable solve($prob)
    g["secant_roots"] = @benchmarkable solve($prob, Secant())
    # Only complex-capable solvers qualify
    # bracketing methods need a real sign change Broyden/DFSane assume a real residual.
    g["muller_sciml"] = @benchmarkable solve(ip, alg) setup = (
        ip=SNS.IntervalNonlinearProblem((ω, _) -> $f(ω), ($prob.omega0 - $h, $prob.omega0 + $h));
        alg=BNS.Muller($prob.omega0 + $h * im))
    g["halley_sciml"] = @benchmarkable solve($prob, SNS.SimpleHalley())
end
