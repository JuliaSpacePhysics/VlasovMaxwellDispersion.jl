# Compile-latency check
# See experiments/QuadGK-seeding/README.md.

@testitem "CoupledVDF cold compile latency vs main" tags = [:latency] begin
    # Measure in a FRESH process: TestItemRunner reuses one worker across testitems,
    proj = dirname(Base.active_project())
    repo = dirname(@__DIR__)        # test/ → repo root
    src = joinpath(repo, "src", "distributions", "CoupledVDF.jl")

    probe = """
    using VlasovMaxwellDispersion
    using VlasovMaxwellDispersion: CoupledVDF, NormalizedSpecies, Wavenumber, contribution
    g0(u, v) = exp(-(u^2 + v^2 + 0.6u * v))
    s = NormalizedSpecies(-1.0, 1.0, CoupledVDF(g0; para = (-8.0, 8.0), perp = 6.0))
    r = @timed contribution(s, 1.2 + 0.05im, Wavenumber(0.3, 0.4))
    print(r.compile_time)
    """
    coldcompile() = parse(Float64, read(`$(Base.julia_cmd()) --project=$proj -e $probe`, String))
    try
        main_src = read(`git -C $repo show main:src/distributions/CoupledVDF.jl`, String)

        ct_cur = coldcompile()
        orig = read(src, String)
        ct_main = try
            write(src, main_src)
            coldcompile()
        finally
            write(src, orig)
        end
        @info "CoupledVDF cold compile_time" current = ct_cur main = ct_main ratio = ct_cur / ct_main
        @test ct_cur < 1.5 * ct_main
    catch
        @test_skip "main:src/distributions/CoupledVDF.jl unavailable (offline/shallow); cannot baseline"
        return
    end

end
