@testitem "NonnegBSpline knot refinement survives residual plateau instead of bisecting data-free cells forever" begin
    μ = 2.0
    f0(u, w) = exp(-μ * sqrt(1 + u^2 + w^2))
    pperp = range(0.0, 5.0, length = 61)
    ppar = range(-5.0, 5.0, length = 121)
    F = [f0(u, w) for w in pperp, u in ppar]
    Fn = F ./ maximum(F)
    fit = fit_grid(NonnegBSpline{3}(rtol = 1.0e-4), pperp, ppar, Fn)
    @test minimum(diff(fit.knots_perp)) > 0
    @test minimum(diff(fit.knots_para)) > 0
    # fit quality survives early termination
    fit_on = [fit(w, u) for w in pperp, u in ppar]
    @test maximum(abs.(fit_on .- Fn)) < 1.0e-3
end
