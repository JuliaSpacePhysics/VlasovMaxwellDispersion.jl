# Validate perp_analytic.jl (cell_bessel_moment/perp_pwpoly) against brute quadrature.
@testitem "perp_analytic moments vs QuadGK" begin
    using QuadGK
    using VlasovMaxwellDispersion: besselj, cell_bessel_moment, perp_pwpoly

    # ∫_vl^vr v * (Σ_k coeffs[k] v^(k-1)) * J_n1(a v) J_n2(a v) dv via quadrature.
    function quad_perp(coeffs, vl, vr, n1, n2, a)
        poly(v) = sum(coeffs[k] * v^(k - 1) for k in eachindex(coeffs))
        f(v) = v * poly(v) * besselj(n1, a * v) * besselj(n2, a * v)
        quadgk(f, vl, vr; rtol = 1e-13)[1]
    end

    cases = [
        # (coeffs, vl, vr, n1, n2, a)
        ([1.0], 0.0, 1.0, 0, 0, 2.0),
        ([0.0, 1.0], 0.0, 1.0, 0, 0, 2.0),
        ([1.0, 2.0, -3.0], 0.5, 1.5, 1, 1, 3.0),
        ([1.0], 0.2, 0.8, 2, 2, 5.0),
        ([1.0, 0.0, 1.0], 0.0, 2.0, 0, 2, 1.5),
        ([1.0, -1.0, 0.5, 0.3], 1.0, 3.0, 3, 1, 4.0),
        ([1.0], 0.0, 5.0, 0, 0, 0.05),                # small a
        ([1.0], 0.0, 1.0, 5, 5, 10.0),                # z=10, near series cutoff
        ([1.0], 0.0, 2.0, 4, 4, 8.0),                 # z=16, BigFloat fallback
        ([1.0], 0.0, 3.0, 0, 0, 15.0),                # z=45, deep BigFloat fallback
        ([2.5, -1.0], 0.1, 0.9, 2, 0, 6.0),
    ]
    for (coeffs, vl, vr, n1, n2, a) in cases
        analytic = cell_bessel_moment(coeffs, vl, vr, n1, n2, a)
        numeric = quad_perp(coeffs, vl, vr, n1, n2, a)
        abserr = abs(analytic - numeric)
        relerr = abserr / max(abs(numeric), 1e-300)
        @test abserr < 1e-8 || relerr < 1e-8
    end

    # perp_pwpoly: multi-cell piecewise polynomial (n1=n2 diagonal moment).
    nodes = [0.0, 0.5, 1.2, 2.0]
    coeffs_cells = [[1.0, 0.5], [0.8, -0.3, 0.1], [0.3]]
    n_, a_ = 2, 3.0
    analytic_pw = perp_pwpoly(coeffs_cells, nodes, n_, n_, a_)
    numeric_pw = sum(quad_perp(coeffs_cells[i], nodes[i], nodes[i + 1], n_, n_, a_) for i in 1:3)
    @test abs(analytic_pw - numeric_pw) / abs(numeric_pw) < 1e-8
end
