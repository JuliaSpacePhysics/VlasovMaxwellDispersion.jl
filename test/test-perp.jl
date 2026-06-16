# Validate perp_analytic.jl (cell_bessel_moment/besselprod_moment) against brute quadrature.
# Not a @testitem: needs QuadGK (test-only dep, not in main Project.toml).
using VlasovMaxwellDispersion
using QuadGK
using Test

const VMD = VlasovMaxwellDispersion
using .VMD: besselj

# ∫_vl^vr v * (Σ_k coeffs[k] v^(k-1)) * J_n1(a v) J_n2(a v) dv via quadrature.
function quad_perp(coeffs, vl, vr, n1, n2, a)
    poly(v) = sum(coeffs[k] * v^(k - 1) for k in eachindex(coeffs))
    f(v) = v * poly(v) * besselj(n1, a * v) * besselj(n2, a * v)
    quadgk(f, vl, vr; rtol=1e-13)[1]
end

cases = [
    # (coeffs, vl, vr, n1, n2, a, label)
    ([1.0], 0.0, 1.0, 0, 0, 2.0, "d=0 n=n=0 a=2 cell[0,1]"),
    ([0.0, 1.0], 0.0, 1.0, 0, 0, 2.0, "d=1 monomial v n=0 a=2"),
    ([1.0, 2.0, -3.0], 0.5, 1.5, 1, 1, 3.0, "quadratic poly n=1 a=3 cell[0.5,1.5]"),
    ([1.0], 0.2, 0.8, 2, 2, 5.0, "d=0 n=2 a=5 cell[0.2,0.8]"),
    ([1.0, 0.0, 1.0], 0.0, 2.0, 0, 2, 1.5, "n1!=n2 (0,2) a=1.5"),
    ([1.0, -1.0, 0.5, 0.3], 1.0, 3.0, 3, 1, 4.0, "cubic poly n1!=n2 (3,1) a=4"),
    ([1.0], 0.0, 5.0, 0, 0, 0.05, "small a (a=0.05) n=0"),
    ([1.0], 0.0, 1.0, 5, 5, 10.0, "high n=5 a=10 -> z=10 (near series cutoff)"),
    ([1.0], 0.0, 2.0, 4, 4, 8.0, "n=4 a=8 -> z=16 (BigFloat fallback path)"),
    ([1.0], 0.0, 3.0, 0, 0, 15.0, "n=0 a=15 -> z=45 (deep BigFloat fallback)"),
    ([2.5, -1.0], 0.1, 0.9, 2, 0, 6.0, "linear poly n1!=n2 (2,0) a=6"),
]

maxerr = 0.0
println("case                                          analytic              quad              abs_err      rel_err")
for (coeffs, vl, vr, n1, n2, a, label) in cases
    analytic = VMD.cell_bessel_moment(coeffs, vl, vr, n1, n2, a)
    numeric = quad_perp(coeffs, vl, vr, n1, n2, a)
    abserr = abs(analytic - numeric)
    relerr = abserr / max(abs(numeric), 1e-300)
    global maxerr = max(maxerr, abserr, relerr)
    println(rpad(label, 45), "  ", analytic, "  ", numeric, "  ", abserr, "  ", relerr)
    @test abserr < 1e-8 || relerr < 1e-8
end

# perp_pwpoly: multi-cell piecewise polynomial (n1=n2 diagonal moment).
nodes = [0.0, 0.5, 1.2, 2.0]
coeffs_cells = [[1.0, 0.5], [0.8, -0.3, 0.1], [0.3]]
n_, a_ = 2, 3.0
analytic_pw = VMD.perp_pwpoly(coeffs_cells, nodes, n_, n_, a_)
numeric_pw = sum(
    quad_perp(coeffs_cells[i], nodes[i], nodes[i+1], n_, n_, a_) for i in 1:3
)
abserr_pw = abs(analytic_pw - numeric_pw)
relerr_pw = abserr_pw / abs(numeric_pw)
global maxerr = max(maxerr, abserr_pw, relerr_pw)
println(rpad("perp_pwpoly multi-cell n=2 a=3", 45), "  ", analytic_pw, "  ", numeric_pw, "  ", abserr_pw, "  ", relerr_pw)
@test abserr_pw < 1e-8 || relerr_pw < 1e-8

println("\nmax(abs_err, rel_err) over all cases = ", maxerr)
