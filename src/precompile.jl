using PrecompileTools: @setup_workload, @compile_workload

# One end-to-end survey caches whole solve stack
# with erased dispersion function for discover, polish, and link layers.
@setup_workload begin
    vdf = BiKappa(vth_para = 4.0e-3, vth_perp = 3.0e-3, kappa = 4)
    plasma = (
        NormalizedSpecies(1.0, 1.0e4, vdf),
        NormalizedSpecies(-1836.15, 1.836e7, Maxwellian(2.0e-2)),
    )
    region = (-0.05 - 0.1im, 0.3 + 0.06im)
    geom = AngleSweep(k = [30.0], theta = deg2rad(45))
    @compile_workload begin
        prob = GlobalDispersionProblem(plasma, region, geom)
        solve(prob, AAA(n = (5, 4)))
    end
end
