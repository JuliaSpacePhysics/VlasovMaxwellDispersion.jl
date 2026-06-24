# Compare different ways to evaluate the σ-quartet (derivation §6.1) that the Newberger kernel `_qin_T_bare` needs:
#
#   series  — `_qin_sigmas_series` : entire z²-series (the package's |z|<1 branch)
#   closed  — `_qin_sigmas_closed` : 4 complex-order Bessel evals (the |z|≥1 branch)
#   pFq     — HypergeometricFunctions: σ0 = (1/a)·₁F₂(½; 1+a, 1−a; −z²)
#
# Reference is the same z²-series in BigFloat (entire ⇒ exact once enough terms are summed
# at high precision).
#
# Run: julia --project=benchmark benchmark/qin_sigmas_compare.jl

# Note: For integer values of a (a ∈ ℤ), the hypergeometric function ₁F₂ has a degenerate parameter.
# HypergeometricFunctions.jl detects it and applies the standard hypergeometric convention
# (regularized / limiting evaluation, Johansson §2.1: "if any bⱼ∈ℤ_{≤0}… conventional to use the truncated series")
# instead of dividing by zero. So it returns a finite number.

using VlasovMaxwellDispersion
using VlasovMaxwellDispersion: _qin_sigmas_series, _qin_sigmas_closed
import HypergeometricFunctions: pFq
using Printf
using Chairmarks: @b

# σ0 via pFq (the others are algebraic consequences; σ0 is the core hypergeometric).
sigma0_pFq(a, z) = (one(a) / a) * pFq((one(real(a)) / 2,), (1 + a, 1 - a), -z^2)

# Gold reference: the σ0 series summed in BigFloat. Entire function ⇒ converges for all z;
# 512-bit precision absorbs the large-z cancellation, so the Float64 truncation is ~1e-16.
function sigma0_ref(a, z)
    return setprecision(BigFloat, 512) do
        ab = Complex{BigFloat}(a)
        zb = BigFloat(z)
        x = (zb / 2)^2
        q = one(ab) / ab
        σ0 = q
        xpow = one(BigFloat)              # x^k starts at k=1 below
        for k in 1:400
            q *= -(2k) * (2k - 1) / (k^2 * (k^2 - ab^2))
            xpow *= x
            term = q * xpow
            σ0 += term
            abs(term) <= eps(BigFloat) * abs(σ0) && k > abs(ab) && break
        end
        ComplexF64(σ0)
    end
end

reld(x, ref) = abs(x - ref) / abs(ref)

function sigma_compare(;
        as = (0.37 + 0.1im, 3.0 + 1.0e-4im),
        zs = (0.1, 0.5, 0.9, 1.0, 2.0, 5.0, 10.0, 30.0)
    )
    for a in as
        @printf("\n=== a = %s   (|a|=%.2f) ===\n", string(a), abs(a))
        println(
            rpad("z", 7), rpad("series err", 13), rpad("closed err", 13), rpad("pFq err", 13),
            rpad("t_series", 11), rpad("t_closed", 11), "t_pFq"
        )
        for z in zs
            ref = sigma0_ref(a, z)
            es = reld(_qin_sigmas_series(a, z)[1], ref)
            ec = reld(_qin_sigmas_closed(a, z)[1], ref)
            ep = reld(sigma0_pFq(a, z), ref)
            ts = (@b _qin_sigmas_series($a, $z)).time
            tc = (@b _qin_sigmas_closed($a, $z)).time
            tp = (@b sigma0_pFq($a, $z)).time
            println(
                rpad(@sprintf("%.2g", z), 7),
                rpad(@sprintf("%.1e", es), 13), rpad(@sprintf("%.1e", ec), 13), rpad(@sprintf("%.1e", ep), 13),
                rpad(@sprintf("%.0f ns", 1.0e9 * ts), 11), rpad(@sprintf("%.0f ns", 1.0e9 * tc), 11),
                @sprintf("%.0f ns", 1.0e9 * tp)
            )
        end
    end
    println("\nerr = rel. error of σ0 vs BigFloat-series gold. Package uses series for |z|<1, closed for |z|≥1.")
    println("Expect: series accurate+fast small z, cancels large z; closed accurate large z, degrades→z=0; pFq uniform but slowest.")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    sigma_compare()
end
