@testitem "survey solve! infers a concrete SurveySolution" begin
    using VlasovMaxwellDispersion: Muller, AAA, CommonSolve

    s = NormalizedSpecies(-1.0, 1.0, Maxwellian(vth_para = 0.9, vth_perp = 1.2))
    prob = GlobalDispersionProblem(s, (0.2 - 0.4im, 1.5 + 0.1im), Wavenumber(0.01, 0.5))
    cache = init(prob, AAA(); refine = Muller())
    rt = only(Base.return_types(CommonSolve.solve!, (typeof(cache),)))
    @test isconcretetype(rt)
end
