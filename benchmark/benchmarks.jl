using BenchmarkTools
using VlasovMaxwellDispersion
using VlasovMaxwellDispersion: contribution
import SimpleNonlinearSolve as SNS
import BracketingNonlinearSolve as BNS
import Roots
const SUITE = BenchmarkGroup()

# Maxwellian with closed Z/Γ_n harmonic sum
let g = SUITE["Maxwellian"]
    s = NormalizedSpecies(-1.0, 1.0, Maxwellian(vth_para=0.9, vth_perp=1.2))
    for kp in (0.1, 1.0)
        k = Wavenumber(kp, 0.4)
        g["kperp=$kp"] = @benchmarkable contribution($s, 1.3 - 0.05im, $k)
    end
end

# Sweep k⊥ — cost grows with nmax ∝ k⊥ρ (harmonic count + per-node Bessel ladder).
let g = SUITE["separable"]
    vthp, vthq = 0.9, 1.2
    sep = SeparableVDF(Maxwellian(vth_para=vthp, vth_perp=vthq); para=(-14vthp, 14vthp), perp=14vthq)
    for kp in (0.1, 1.0, 2.5)
        k = Wavenumber(kp, 0.4)
        g["gaussian_kperp=$kp"] = @benchmarkable contribution($sep, 1.3 - 0.05im, $k)
    end
    # Non-Gaussian parallel (kappa-like) × Gaussian perp
    sk = SeparableVDF(v -> exp(-v^2) / pi, u -> (1 + u^2 / 3)^(-2); para=(-30.0, 30.0), perp=10.0)
    for kp in (1.0, 3.0)
        k = Wavenumber(kp, 0.4)
        g["kappa_kperp=$kp"] = @benchmarkable contribution($sk, 1.2 - 0.05im, $k)
    end
end

let g = SUITE["kappa"]
    ω = 1.2 - 0.05im
    for κ in (6, 6.00001)   # integer (residue) vs non-integer (₂F₁) parallel path
        s = NormalizedSpecies(-1.0, 0.7, BiKappa(vth_para=0.9, vth_perp=1.2, kappa=κ))
        for kp in (0.4, 2.0)
            g["bikappa=$κ/kperp=$kp"] = @benchmarkable contribution($s, $ω, Wavenumber($kp, 0.3))
        end
    end
    pbk = NormalizedSpecies(-1.0, 0.7, ProductBiKappa(vth_para=0.9, vth_perp=1.2, kappa_para=6, kappa_perp=4))
    for kp in (0.4, 2.0)
        g["product/kperp=$kp"] = @benchmarkable contribution($pbk, $ω, Wavenumber($kp, 0.3))
    end
end

# CoupledVDF: the two closures (derivation.md §3) on an inseparable f₀.
let g = SUITE["Nonrelativistic/coupled"]
    g0(w, u) = exp(-(u^2 + w^2 + 0.6u * w))
    kw = (para=(-8.0, 8.0), perp=6.0)
    s = NormalizedSpecies(-1.0, 1.0, CoupledVDF(g0; kw...))
    ω = 1.2 + 0.05im
    for kp in (0.3, 1.0)
        k = Wavenumber(kp, 0.4)
        g["B_truncated/kperp=$kp"] = @benchmarkable contribution($s, $ω, $k)
        g["A_newberger/kperp=$kp"] = @benchmarkable contribution($s, $ω, $k; closure=Newberger())
    end
end

let g = SUITE["Relativistic"], regime=Relativistic()
    μ = 40.0
    γ(w, u) = sqrt(1 + u^2 + w^2)
    f0(w, u) = exp(-μ * γ(w, u))
    L = sqrt((1 + 25 / μ)^2 - 1)
    kw = (para=(-L, L), perp=L)
    vdf = CoupledVDF(f0; kw..., regime)
    ω = 0.3 + 0.05im
    # Sweep k⊥: edge-mapped (γ,p∥) cost grows with nmax∝k⊥ρ
    for kp in (0.7, 3.5)
        k = Wavenumber(kp, 0.4)
        g["coupled/A_newberger/kperp=$kp"] = @benchmarkable contribution($vdf, $ω, $k; closure=Newberger())
        g["coupled/B_truncated/kperp=$kp"] = @benchmarkable contribution($vdf, $ω, $k)
    end
    k = Wavenumber(0.7, 0.4)

    ppar = collect(range(-L, L, length=81))
    pperp = collect(range(0.0, L, length=61))
    F = [f0(u, w) for w in pperp, u in ppar]      # F[perp,par]
    vdf = GridVDF(pperp, ppar, F; regime)
    g["gridvdf"] = @benchmarkable contribution($vdf, $ω, $k)
    g["gridvdf_coupled"] = @benchmarkable contribution($vdf.coupled, $ω, $k)
end

let g = SUITE["local_solve"]
    s = NormalizedSpecies(-1.0, 1.0, Maxwellian(vth_para=0.9, vth_perp=1.2))
    k = Wavenumber(0.01, 0.5)
    prob = LocalDispersionProblem(s, k, 0.6)
    h = 1e-3
    g["muller_native"] = @benchmarkable solve($prob)
    g["secant_roots"] = @benchmarkable solve($prob, Roots.Order1())
    # Only complex-capable solvers qualify
    # bracketing methods need a real sign change Broyden/DFSane assume a real residual.
    g["muller_sciml"] = @benchmarkable solve(ip, alg) setup = (
        ip=SNS.IntervalNonlinearProblem((ω, _) -> $prob.f(ω), ($prob.omega0 - $h, $prob.omega0 + $h));
        alg=BNS.Muller($prob.omega0 + $h * im))
    g["halley_sciml"] = @benchmarkable solve($prob, SNS.SimpleHalley())
end
