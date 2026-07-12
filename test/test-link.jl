using VlasovMaxwellDispersion: link

@testset "link joins one constant root into a single m=2 sheet" begin
    vals = fill([1.0 + 0.0im], 5, 5)
    sheets = link(vals; gate = 0.5)
    @test length(sheets) == 1
    @test all(==(1.0 + 0.0im), only(sheets))
end

@testset "link separates far values into distinct sheets" begin
    vals = [[1.0 + 0.0im, 5.0 + 0.0im] for _ in 1:5]
    sheets = link(vals; gate = 0.5)
    @test length(sheets) == 2
    @test all(sh -> count(isfinite, sh) == 5, sheets)
end
