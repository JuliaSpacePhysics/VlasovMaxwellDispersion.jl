# χ(ω) is entire for an entire f₀, so ω²χ obeys the mean-value property on any circle. 
# Fixed-node plans switch branches inside such circles
# (residue on/off at the real axis, pole peeling near it); an inconsistent switch breaks the identity.
# Circles: across the real axis (peel and residue switch on), and deep
# on the damped side at small k∥, where crossed poles lie beyond the peel band.
@testitem "Fixed-node plans: ω²χ is holomorphic across their branch switches" begin
    using VlasovMaxwellDispersion: plan_contribution
    g0(w, u) = exp(-(u^2 + w^2 + 0.6u * w))
    vdfs = (
        CoupledVDF(g0; para = 8.0, perp = 6.0),
        SeparableVDF(v -> exp(-v^2), u -> (1 + u^2 / 2) * exp(-u^2); para = 8.0, perp = 6.0),
        LowRankVDF(g0; para = 8.0, perp = 6.0, rtol = 1.0e-10),
    )
    meanvalue(f, ω0, r; N = 256) = sum(j -> f(ω0 + r * cis(2π * j / N)), 0:(N - 1)) / N
    circles = (
        (Wavevector(4.0, 0.4), 1.3 + 0.0im, 0.1),
        (Wavevector(4.0, 0.4), 0.6 + 0.02im, 0.5),
        (Wavevector(0.3, -0.4), 1.3 - 0.05im, 0.3),
        (Wavevector(1.0, 0.25), 1.0 - 0.6im, 0.5),
    )
    for vdf in vdfs
        s = prepare(NormalizedPopulation(-1.0, 1.0, vdf))
        for (k, ω0, r) in circles
            pl = plan_contribution(s, k)
            f(ω) = ω^2 * pl(ω)
            @test meanvalue(f, ω0, r) ≈ f(ω0) rtol = 1.0e-6
        end
    end
end
