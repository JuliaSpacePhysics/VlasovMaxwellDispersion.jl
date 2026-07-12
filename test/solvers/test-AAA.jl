@testitem "AAA finds near-origin roots the raw det's ω=0 pole would screen" begin
    using VlasovMaxwellDispersion: NormalizedSpecies, Wavenumber, ColdVDF,
        GlobalDispersionProblem, AAA, solve

    # Cold e-p plasma whose parallel modes sit near ω=0 (dense: low vA).
    mp_me = 1836.15
    plasma = (
        NormalizedSpecies(-1.0, 100.0, ColdVDF()),
        NormalizedSpecies(1 / mp_me, 100.0 / mp_me, ColdVDF()),
    )
    k = Wavenumber(0.0, 2.0)
    region = (-0.08 - 0.03im, 0.08 + 0.03im)   # straddles ω=0
    prob = GlobalDispersionProblem(plasma, region, k)

    sol = solve(prob, AAA())
    @test sol.retcode == :Success

    # The genuine symmetric ±kvA root pair is found off the origin
    @test count(b -> abs(b.omega) > 0.02, sol.roots) ≥ 2
    @test all(b.resid < 1.0e-4 for b in sol.roots)
end

@testitem "AAA per-slice (m=1) tracks a root branch" begin
    using VlasovMaxwellDispersion: NormalizedSpecies, ColdVDF, CartesianSweep,
        GlobalDispersionProblem, AAA, GRPF, solve, para

    cold = (
        NormalizedSpecies(-1.0, 10.0, ColdVDF()),
        NormalizedSpecies(1 / 1836.15, 10.0 / 1836.15, ColdVDF()),
    )
    region = (0.05 - 0.05im, 1.2 + 0.05im)
    prob = GlobalDispersionProblem(cold, region, CartesianSweep(; kz = (0.3, 3.0)))

    sp = solve(prob, AAA())
    @test sp.retcode == :Success
    @test !isempty(sp.roots)
    rb = argmax(b -> count(isfinite, b.omega), sp.roots)
    valid = isfinite.(rb.omega)
    # Every sample Muller-polished to a genuine zero of det𝒟.
    @test maximum(rb.resid[valid]) < 1.0e-6
    # Linked into one sheet spanning the bulk of the sweep, not per-k fragments.
    kz = para.(rb.k)
    @test maximum(kz[valid]) - minimum(kz[valid]) > 1.5

    # Cross-check a mid-sweep sample against an independent fixed-k GRPF survey.
    i = argmin(i -> valid[i] ? abs(kz[i] - 1.5) : Inf, eachindex(kz))
    sg = solve(GlobalDispersionProblem(cold, region, rb.k[i]), GRPF())
    @test minimum(abs(b.omega - rb.omega[i]) for b in sg.roots) < 1.0e-6
end

@testitem "AAA and padded GRPF swept surveys agree" begin
    using VlasovMaxwellDispersion: NormalizedSpecies, ColdVDF, CartesianSweep,
        GlobalDispersionProblem, AAA, GRPF, solve, para

    cold = (
        NormalizedSpecies(-1.0, 10.0, ColdVDF()),
        NormalizedSpecies(1 / 1836.15, 10.0 / 1836.15, ColdVDF()),
    )
    # Box holds four parallel branches: whistler, two EM branches, and the
    # k-independent Langmuir mode at ω = ωpe = √10.005 ≈ 3.163.
    # The σ_min/σ_max residual stays ~ε (asserted below).
    region = (0.05 - 0.05im, 4.5 + 0.05im)
    prob = GlobalDispersionProblem(cold, region, CartesianSweep(; kz = range(0.3, 3.0, length = 9)))
    sa = solve(prob, AAA())
    sg_exact = solve(prob, GRPF())
    ll, ur = region
    margin = 0.15 * (ur - ll)
    padded = (ll - margin, ur + margin)
    sg_padded = solve(GlobalDispersionProblem(cold, padded, prob.geometry), GRPF())
    @test sa.retcode == sg_exact.retcode == sg_padded.retcode == :Success
    @test length(sa.roots) == 4

    inbox(ω) = real(ll) ≤ real(ω) ≤ real(ur) && imag(ll) ≤ imag(ω) ≤ imag(ur)
    samples(sol) = [(para(k), ω) for b in sol.roots for (k, ω) in zip(b.k, b.omega) if inbox(ω)]
    sA, sG_exact, sG = samples(sa), samples(sg_exact), samples(sg_padded)
    # GRPF misses two roots near boundaries
    @test length(sG_exact) < length(sA) == length(sG)

    dist(x, S) = minimum(hypot(x[1] - y[1], abs(x[2] - y[2])) for y in S)
    @test maximum(dist(x, sG) for x in sA) < 1.0e-8
    @test maximum(dist(x, sA) for x in sG) < 1.0e-8
    @test any(b -> all(ω -> abs(ω - sqrt(10.005446)) < 1.0e-3, filter(isfinite, b.omega)), sa.roots)
    # σ-ratio residual is ~ε on every branch, including row-vanishing Langmuir.
    @test all(maximum(filter(isfinite, b.resid)) < 1.0e-6 for b in sa.roots)
end

@testitem "AAA keeps genuine low-frequency sheets in a wide ω box" begin
    using VlasovMaxwellDispersion: NormalizedSpecies, ColdVDF, AngleSweep,
        GlobalDispersionProblem, AAA, SurveySolution, solve

    # Electron-scale ω box over an e-p plasma: the ion-cyclotron branch lives at
    # ω ≲ ωci ≈ 5.4e-4·ωce — three decades below the box scale. A magnitude-based
    # artifact heuristic would delete the whole sheet.
    mp_me = 1836.15
    plasma = (
        NormalizedSpecies(-1.0, 100.0, ColdVDF()),
        NormalizedSpecies(1 / mp_me, 100.0 / mp_me, ColdVDF()),
    )
    region = (1.0e-4 - 0.05im, 3.0 + 0.01im)
    # Sweep resolution is the geometry's: an explicit k grid replaces any solver knob.
    geom = AngleSweep(k = range(0.02, 2.0, length = 31), theta = 0.001)
    sol = solve(GlobalDispersionProblem(plasma, region, geom), AAA())
    @test sol.retcode == :Success
    @test any(b -> maximum(abs, filter(isfinite, b.omega)) < 3.0e-3, sol.roots)
    @test all(maximum(filter(isfinite, b.resid)) < 1.0e-4 for b in sol.roots)

    # Post-hoc pruning without re-solving.
    long = filter(b -> count(isfinite, b.omega) ≥ 10, sol)
    @test long isa SurveySolution
    @test all(count(isfinite, b.omega) ≥ 10 for b in long.roots)
end

@testitem "AAA rejects the deflated det's structural ω=0 zero (kinetic)" begin
    using VlasovMaxwellDispersion: NormalizedSpecies, Wavenumber, ProductBiKappa,
        GlobalDispersionProblem, AAA, solve, residual, wave_dispersion_tensor
    using LinearAlgebra: det

    # Kinetic species ⇒ ω²χ → 0 at ω=0 ⇒ det(ω²𝒟) has a zero pinned at the
    # origin. The raw residual CANNOT reject it (it also vanishes as ω → 0) — only
    # the geometric |ω| gate can.
    vdf = ProductBiKappa(vth_para = 0.1, vth_perp = 0.1, kappa_para = 1, kappa_perp = 200.0)
    plasma = (NormalizedSpecies(-1.0, 300.0, vdf),)
    k = Wavenumber(0.5, 1.0)
    @test abs(det(wave_dispersion_tensor(plasma, 1.0e-5 + 0im, k))) < 1.0e-4
    @test residual(plasma, 1.0e-5 + 0im, k) < 1.0e-4   # residual gate blind here

    region = (-0.02 - 0.02im, 0.35 + 0.01im)           # straddles the origin
    sol = solve(GlobalDispersionProblem(plasma, region, k), AAA())
    diag = hypot(0.37, 0.03)
    @test all(abs(b.omega) > 1.0e-6 * diag for b in sol.roots)
    @test all(b.resid < 1.0e-4 for b in sol.roots)
end
