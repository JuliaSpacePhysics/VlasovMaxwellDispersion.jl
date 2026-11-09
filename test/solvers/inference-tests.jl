@testitem "survey solve! infers a concrete SurveySolution" begin
    using VlasovMaxwellDispersion: init, solve!, SurveySolution
    using VlasovMaxwellDispersion: Muller, AAA, CommonSolve

    s = NormalizedPopulation(-1.0, 1.0, BiMaxwellian(; vth_perp = 1.2, vth_para = 0.9))
    prob = DispersionProblem(s, (0.2 - 0.4im, 1.5 + 0.1im), Wavevector(0.01, 0.5))
    cache = init(prob, AAA(); refine = Muller())
    rt = only(Base.return_types(CommonSolve.solve!, (typeof(cache),)))
    @test isconcretetype(rt)
end
