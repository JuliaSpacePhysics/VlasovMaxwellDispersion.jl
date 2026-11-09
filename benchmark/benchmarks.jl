using ChairmarksForAirspeedVelocity
import ChairmarksForAirspeedVelocity as CAV
using VlasovMaxwellDispersion
using VlasovMaxwellDispersion: contribution, plan_contribution
using PlasmaPopulations: pdf
const SUITE = BenchmarkGroup()

# A solve costs (χ evaluations) × (cost per evaluation), and a plan change can move the two in
# opposite directions, so every plan-backed VDF carries a `global_solve` beside its `evaluate`.
const REGION = (0.2 - 0.4im, 1.5 + 0.1im)

# Maxwellian with closed Z/Γ_n harmonic sum
let g = SUITE["Maxwellian"]
    s = NormalizedPopulation(-1.0, 1.0, BiMaxwellian(; vth_perp = 1.2, vth_para = 0.9))
    g["kperp=1.0"] = CAV.@benchmarkable contribution($s, 1.3 - 0.05im, $(Wavevector(1.0, 0.4)))
end

let g = SUITE["separable"]
    vthp, vthq = 0.9, 1.2
    ssep = prepare(NormalizedPopulation(-1.0, 1.0, SeparableVDF(
        Gaussian(vthq) ⊗ Gaussian(vthp, 0.0);
        para=14vthp, perp=14vthq
    )))
    k = Wavevector(4.0, 0.4)
    g["plan"] = CAV.@benchmarkable plan_contribution($ssep, $k)
    plan = plan_contribution(ssep, k)
    g["evaluate"] = CAV.@benchmarkable $plan(1.3 - 0.05im)
    prob = DispersionProblem(ssep, REGION, Wavevector(1.0, 0.4))
    g["global_solve"] = CAV.@benchmarkable solve($prob, AAA())
end

let g = SUITE["kappa"]
    ω = 1.2 - 0.05im
    k = Wavevector(2.0, 0.3)
    for κ in (6, 6.00001)   # integer (residue) vs non-integer (₂F₁) parallel path
        s = NormalizedPopulation(-1.0, 0.7, BiKappa(; vth_perp = 1.2, vth_para = 0.9, kappa = κ))
        g["bikappa=$κ"] = CAV.@benchmarkable contribution($s, $ω, $k)
    end
    pbk = prepare(NormalizedPopulation(-1.0, 0.7,
        ProductBiKappa(; vth_perp = 1.2, vth_para = 0.9, kappa_perp = 4, kappa_para = 6)
    ))
    g["product/plan"] = CAV.@benchmarkable plan_contribution($pbk, $k)
    plan = plan_contribution(pbk, k)
    g["product/evaluate"] = CAV.@benchmarkable $plan($ω)
end

# CoupledVDF on an inseparable f₀: the default fixed-node plan, and the per-ω adaptive path
# under both closures (derivation.md §3).
let g = SUITE["Nonrelativistic/coupled"]
    g0(w, u) = exp(-(u^2 + w^2 + 0.6u * w))
    s = NormalizedPopulation(-1.0, 1.0, CoupledVDF(g0; para=8.0, perp=6.0))
    g["prepare"] = CAV.@benchmarkable prepare($s)
    sp = prepare(s)
    k = Wavevector(4.0, 0.4)
    g["plan"] = CAV.@benchmarkable plan_contribution($sp, $k)
    plan = plan_contribution(sp, k)
    g["evaluate/undamped"] = CAV.@benchmarkable $plan(1.2 + 0.05im)
    g["evaluate/damped"] = CAV.@benchmarkable $plan(1.2 - 0.05im)   # adds the exact Landau residue
    prob = DispersionProblem(sp, REGION, Wavevector(1.0, 0.4))
    g["global_solve"] = CAV.@benchmarkable solve($prob, AAA())
    k1 = Wavevector(1.0, 0.4)
    g["adaptive/B_truncated"] = CAV.@benchmarkable contribution($sp, 1.2 + 0.05im, $k1)
    g["adaptive/A_newberger"] = CAV.@benchmarkable contribution($sp, 1.2 + 0.05im, $k1; closure=Newberger())
end

let g = SUITE["Nonrelativistic/lowrank"]
    vthz, vthp = 0.2, 0.3
    bk = BiKappa(; vth_perp = vthp, vth_para = vthz, kappa = 3.0)
    kw = (para=20vthz, perp=20vthp, rtol=1.0e-10)
    g["construct"] = CAV.@benchmarkable LowRankVDF($((q, u) -> pdf(bk, q, u)); $kw...)

    vdf = LowRankVDF((q, u) -> pdf(bk, q, u); kw...)
    s = NormalizedPopulation(1.0, 1.0, vdf)
    k = Wavevector(8.0, 0.1)
    g["plan"] = CAV.@benchmarkable plan_contribution($s, $k)
    plan = plan_contribution(s, k)
    g["evaluate"] = CAV.@benchmarkable $plan(1.3 - 0.02im)
    prob = DispersionProblem(s, REGION, Wavevector(2.0, 1.0))
    g["global_solve"] = CAV.@benchmarkable solve($prob, AAA())
end

let g = SUITE["Relativistic"], regime = Relativistic()
    γ(w, u) = sqrt(1 + u^2 + w^2)
    kp = 0.7
    for μ in (2.0, 40.0)
        f0(w, u) = exp(-μ * γ(w, u))
        L = sqrt((1 + 25 / μ)^2 - 1)
        vdf = prepare(CoupledVDF(f0; para=L, perp=L, regime))
        ω = 0.3 + 0.05im
        k = Wavevector(kp, 0.4)
        g["coupled/A_newberger/μ=$μ/kperp=$kp"] = CAV.@benchmarkable contribution($vdf, $ω, $k; closure=Newberger())
        g["coupled/B_truncated/μ=$μ/kperp=$kp"] = CAV.@benchmarkable contribution($vdf, $ω, $k)

        ppar = collect(range(-L, L, length=81))
        pperp = collect(range(0.0, L, length=61))
        F = [f0(u, w) for w in pperp, u in ppar]      # F[perp,par]
        vdf = SplineVDF(pperp, ppar, F; regime)
        g["gridvdf/μ=$μ/kperp=$kp"] = CAV.@benchmarkable contribution($vdf, $ω, $k)
    end
end

let g = SUITE["local_solve"]
    s = NormalizedPopulation(-1.0, 1.0, BiMaxwellian(; vth_perp = 1.2, vth_para = 0.9))
    k = Wavevector(0.01, 0.5)
    prob = DispersionProblem(s, 0.6, k)
    g["muller_native"] = CAV.@benchmarkable solve($prob, Muller())
end

let g = SUITE["branch_solve"]
    s = NormalizedPopulation(-1.0, 1.0, BiMaxwellian(; vth_perp = 1.2, vth_para = 0.9))
    ks = [Wavevector(0.01, kz) for kz in range(0.3, 0.8; length=64)]
    prob = DispersionProblem(s, 0.6, ks)
    g["arc_length/64"] = CAV.@benchmarkable solve($prob, Continuation())
end

let g = SUITE["global_solve"]
    s = NormalizedPopulation(-1.0, 1.0, BiMaxwellian(; vth_perp = 1.2, vth_para = 0.9))
    prob = DispersionProblem(s, REGION, Wavevector(0.01, 0.5))
    g["Default/fixed_k"] = CAV.@benchmarkable solve($prob, AAA())
end
