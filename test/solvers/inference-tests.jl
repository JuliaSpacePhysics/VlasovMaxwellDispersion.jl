@testitem "solve infers concrete solutions" begin
    s = NormalizedPopulation(-1.0, 1.0, BiMaxwellian(; vth_perp = 1.2, vth_para = 0.9))
    k = Wavevector(0.01, 0.5)
    infers(prob, alg) = isconcretetype(only(Base.return_types(solve, (typeof(prob), typeof(alg)))))
    @test infers(DispersionProblem(s, (0.2 - 0.4im, 1.5 + 0.1im), k), AAA())
    @test infers(DispersionProblem(s, 0.6, k), Muller())
    @test infers(DispersionProblem(s, 0.6, [k]), Continuation())
end
