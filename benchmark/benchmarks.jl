using ChairmarksForAirspeedVelocity
import ChairmarksForAirspeedVelocity as CAV
using VlasovMaxwellDispersion
using VlasovMaxwellDispersion: contribution
import SimpleNonlinearSolve as SNS
import BracketingNonlinearSolve as BNS
const SUITE = BenchmarkGroup()

# Maxwellian with closed Z/Γ_n harmonic sum
let g = SUITE["Maxwellian"]
    s = NormalizedSpecies(-1.0, 1.0, Maxwellian(vth_para = 0.9, vth_perp = 1.2))
    for kp in (0.1, 1.0)
        k = Wavenumber(kp, 0.4)
        g["kperp=$kp"] = CAV.@benchmarkable contribution($s, 1.3 - 0.05im, $k)
    end
end

# Sweep k⊥ — cost grows with nmax ∝ k⊥ρ (harmonic count + per-node Bessel ladder).
let g = SUITE["separable"]
    vthp, vthq = 0.9, 1.2
    sep = prepare(SeparableVDF(Maxwellian(vth_para = vthp, vth_perp = vthq); para = (-14vthp, 14vthp), perp = 14vthq))
    for kp in (0.1, 1.0, 2.5)
        k = Wavenumber(kp, 0.4)
        g["gaussian_kperp=$kp"] = CAV.@benchmarkable contribution($sep, 1.3 - 0.05im, $k)
    end
    # Non-Gaussian parallel (kappa-like) × Gaussian perp
    sk = prepare(SeparableVDF(v -> exp(-v^2) / pi, u -> (1 + u^2 / 3)^(-2); para = (-30.0, 30.0), perp = 10.0))
    for kp in (1.0, 3.0)
        k = Wavenumber(kp, 0.4)
        g["kappa_kperp=$kp"] = CAV.@benchmarkable contribution($sk, 1.2 - 0.05im, $k)
    end
    # Fixed-k plan lifecycle: hoisting the perp Bessel tensors out of the ω loop.
    ssep = NormalizedSpecies(-1.0, 1.0, SeparableVDF(Maxwellian(vth_para = vthp, vth_perp = vthq); para = (-14vthp, 14vthp), perp = 14vthq))
    for (name, k) in (("moderate", Wavenumber(1.0, 0.4)), ("many_harmonics", Wavenumber(4.0, 0.1)))
        g["plan/$name"] = CAV.@benchmarkable plan_contribution($ssep, $k)
        plan = plan_contribution(ssep, k)
        g["evaluate/$name"] = CAV.@benchmarkable $plan(1.3 - 0.05im)
    end
end

let g = SUITE["kappa"]
    ω = 1.2 - 0.05im
    for κ in (6, 6.00001)   # integer (residue) vs non-integer (₂F₁) parallel path
        s = NormalizedSpecies(-1.0, 0.7, BiKappa(vth_para = 0.9, vth_perp = 1.2, kappa = κ))
        for kp in (0.4, 2.0)
            g["bikappa=$κ/kperp=$kp"] = CAV.@benchmarkable contribution($s, $ω, Wavenumber($kp, 0.3))
        end
    end
    pbk = NormalizedSpecies(-1.0, 0.7, ProductBiKappa(vth_para = 0.9, vth_perp = 1.2, kappa_para = 6, kappa_perp = 4))
    for kp in (0.4, 2.0)
        g["product/kperp=$kp"] = CAV.@benchmarkable contribution($pbk, $ω, Wavenumber($kp, 0.3))
    end
end

# CoupledVDF: the two closures (derivation.md §3) on an inseparable f₀.
let g = SUITE["Nonrelativistic/coupled"]
    g0(w, u) = exp(-(u^2 + w^2 + 0.6u * w))
    kw = (para = (-8.0, 8.0), perp = 6.0)
    s = NormalizedSpecies(-1.0, 1.0, CoupledVDF(g0; kw...))
    sp = prepare(s)
    ω = 1.2 + 0.05im
    for kp in (0.3, 1.0, 3.0)
        k = Wavenumber(kp, 0.4)
        g["B_truncated/kperp=$kp"] = CAV.@benchmarkable contribution($sp, $ω, $k)
        g["A_newberger/kperp=$kp"] = CAV.@benchmarkable contribution($sp, $ω, $k; closure = Newberger())
    end
end

let g = SUITE["Nonrelativistic/lowrank"]
    vthz, vthp = 0.2, 0.3
    bk = BiKappa(vth_para = vthz, vth_perp = vthp, kappa = 3.0)
    kw = (para = (-20vthz, 20vthz), perp = 20vthp, rtol = 1.0e-10)
    g["construct"] = CAV.@benchmarkable LowRankVDF($bk; $kw...)

    vdf = LowRankVDF(bk; kw...)
    s = NormalizedSpecies(1.0, 1.0, vdf)
    for (name, k) in (("moderate", Wavenumber(2.0, 1.0)), ("many_harmonics", Wavenumber(8.0, 0.1)))
        g["plan/$name"] = CAV.@benchmarkable plan_contribution($s, $k)
        plan = plan_contribution(s, k)
        g["evaluate/$name"] = CAV.@benchmarkable $plan(1.3 - 0.02im)
    end
    prob = DispersionProblem(s, (0.2 - 0.4im, 1.5 + 0.1im), Wavenumber(2.0, 1.0))
    g["global_solve/moderate"] = CAV.@benchmarkable solve($prob)
end

let g = SUITE["Relativistic"], regime = Relativistic()
    γ(w, u) = sqrt(1 + u^2 + w^2)
    for μ in (2.0, 40.0), kp in (0.7, 3.5)
        f0(w, u) = exp(-μ * γ(w, u))
        L = sqrt((1 + 25 / μ)^2 - 1)
        vdf = prepare(CoupledVDF(f0; para = (-L, L), perp = L, regime))
        ω = 0.3 + 0.05im
        k = Wavenumber(kp, 0.4)
        g["coupled/A_newberger/μ=$μ/kperp=$kp"] = CAV.@benchmarkable contribution($vdf, $ω, $k; closure = Newberger())
        g["coupled/B_truncated/μ=$μ/kperp=$kp"] = CAV.@benchmarkable contribution($vdf, $ω, $k)

        ppar = collect(range(-L, L, length = 81))
        pperp = collect(range(0.0, L, length = 61))
        F = [f0(u, w) for w in pperp, u in ppar]      # F[perp,par]
        vdf = GridVDF(pperp, ppar, F; regime)
        g["gridvdf/μ=$μ/kperp=$kp"] = CAV.@benchmarkable contribution($vdf, $ω, $k)
        g["gridvdf_coupled/μ=$μ/kperp=$kp"] = CAV.@benchmarkable contribution($vdf.coupled, $ω, $k)
    end
end

let g = SUITE["local_solve"]
    s = NormalizedSpecies(-1.0, 1.0, Maxwellian(vth_para = 0.9, vth_perp = 1.2))
    k = Wavenumber(0.01, 0.5)
    prob = DispersionProblem(s, 0.6, k)
    h = 1.0e-3
    g["muller_native"] = CAV.@benchmarkable solve($prob)
    # Only complex-capable solvers qualify
    # bracketing methods need a real sign change Broyden/DFSane assume a real residual.
    g["muller_sciml"] = CAV.@benchmarkable (
        ip = SNS.IntervalNonlinearProblem((ω, _) -> $prob.f(ω), ($prob.omega0 - $h, $prob.omega0 + $h));
        alg = BNS.Muller($prob.omega0 + $h * im);
        (ip, alg)
    ) (x -> solve(x...))
    g["halley_sciml"] = CAV.@benchmarkable solve($prob, SNS.SimpleHalley())
end


let g = SUITE["branch_solve"]
    s = NormalizedSpecies(-1.0, 1.0, Maxwellian(vth_para = 0.9, vth_perp = 1.2))
    ks = [Wavenumber(0.01, kz) for kz in range(0.3, 0.8; length = 64)]
    prob = DispersionProblem(s, 0.6, ks)
    g["arc_length/64"] = CAV.@benchmarkable solve($prob)
end


let g = SUITE["global_solve"]
    s = NormalizedSpecies(-1.0, 1.0, Maxwellian(vth_para = 0.9, vth_perp = 1.2))
    prob = DispersionProblem(s, (0.2 - 0.4im, 1.5 + 0.1im), Wavenumber(0.01, 0.5))
    g["Default/fixed_k"] = CAV.@benchmarkable solve($prob)
end
