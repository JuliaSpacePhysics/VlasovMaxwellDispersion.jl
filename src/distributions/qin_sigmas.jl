const _QIN_ZC2 = 64.0
@inline qin_sigmas(a, z) =
    abs2(z) < _QIN_ZC2 ? _qin_sigmas_series(a, z) : _qin_sigmas_closed(a, z)

function _qin_sigmas_closed(a, z)
    Ja, Jad = besselj_deriv(a, z)
    J_a, J_ad = besselj_deriv(-a, z)
    s = sinpi(a)
    z2 = z^2
    σ0 = π * J_a * Ja / s
    σ1 = (a * σ0 - one(σ0)) / z2
    σD = ((z / 2) * π * (J_ad * Ja + J_a * Jad) / s) / z2
    σJ = π * J_ad * Jad / s + a / z2
    return σ0, σ1, σD, σJ
end

# Series branch: the entire z²-series with q_k≡(π/sin πa)·p_k and the product-series recurrence
# p_k/p_{k-1} = −(2k)(2k−1)/(k²(k²−a²)) ⇒ q_k = q_{k-1}·−(2k)(2k−1)/(k²(k²−a²))
# Two convergence guards, both from Johansson, *Computing Hypergeometric Functions Rigorously* (Thm 1 tail bound):
# (i) gate `k>|a|` so the stop never fires on the pre-convergent terms before the `k≈|a|` spike of the `1/(k²−a²)` resonance;
# (ii) test the k²-weighted σJ increment (the slowest of the four sums) against the `1/|a|` scale of σ0 —
#      once past the spike the term ratio is <1, so a negligible k²-term bounds the geometric tail of every accumulator.
function _qin_sigmas_series(a, z)
    x = (z / 2)^2                                   # P is a series in (z/2)²
    absa = abs(a)
    tolscale = eps(real(z)) / absa                  # tol · |σ0 leading| = tol/|a|
    q = one(a) / a                                  # q_0 = 1/a
    σ0 = q                                          # Σ q_k x^k
    σ1acc = zero(a)                                 # Σ_{k≥1} q_k x^{k-1}
    σDacc = zero(a)                                 # Σ_{k≥1} k q_k x^{k-1}
    σJacc = zero(a)                                 # Σ_{k≥1} k² q_k x^{k-1}
    xpow = one(real(z))                             # x^{k-1}
    k = 1
    converged = false
    while k <= 100
        q *= -(2k) * (2k - 1) / (k^2 * (k^2 - a^2))
        term = q * xpow                             # q_k x^{k-1}
        σ1acc += term
        σDacc += k * term
        σJacc += k^2 * term
        σ0 += q * (xpow * x)                        # q_k x^k
        if k > absa && k^2 * abs(term * x) <= tolscale
            converged = true
            break
        end
        xpow *= x
        k += 1
    end
    converged || return _qin_sigmas_closed(a, z)
    σ1 = (a / 4) * σ1acc
    σD = σDacc / 4
    σJ = σJacc / 2 + σ0 - a * σ1                     # (½Σk²q x^{k-1}) + σ0 − a σ1
    return σ0, σ1, σD, σJ
end
