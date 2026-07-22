# Transported residue cycles: subluminal-germ continuation of the relativistic
# harmonic susceptibility at damped ω (any k⊥, k∥ ≠ 0), generic over the VDF
#   χ_germ = χ_straight − Σₙ (2πi/|k∥|) ∫_Γₙ 𝒰(γ) 𝓣ₙ(aqₙ,uₙ,qₙ) dγ,
#   𝒰 = ω ∂_γf + k∥ ∂_∥f on shell,  uₙ = (γω−nΩ)/k∥,  qₙ² = γ²−1−uₙ².
# 𝓣ₙ is even+entire in qₙ, so for entire denergy the cycle integrand is ENTIRE
# in γ: Γₙ is pinned by (endpoint γ₀ₙ, tail sector) alone.
# Math + certification: docs/src/relativistic.typ.

# cycle endpoint: entering root of qₙ²(γ) = 0; principal Dₙ = √(N²+k∥²−ω²) IS
# the subluminal-germ branch.
@inline function _cycle_endpoint_D(N, ω, kz)
    akz = abs(kz)
    D = sqrt(complex(N^2 + akz^2 - ω^2))
    den = N * ω + akz * D
    γ0 = abs(den) > 1.0e-3 * (abs(N * ω) + abs(akz * D)) ? (N^2 + akz^2) / den :
        (N * ω - akz * D) / (ω^2 - akz^2)   # forms degenerate together only on the light line (gated)
    return γ0, D
end

@inline _cycle_endpoint(N, ω, kz) = first(_cycle_endpoint_D(N, ω, kz))

# Use the known endpoint and path increment: forming γ²−1−u² loses both near the light line.
@inline function _cycle_q2(γ, γ0, u0, ω, kz; D = nothing, δ = nothing)
    δ = @something δ γ - γ0
    slope = isnothing(D) ? 2 * (γ0 - (ω / kz) * u0) : 2D / abs(kz)
    return δ * (slope + (1 - (ω / kz)^2) * δ)
end

# tail direction: maximize e^{−μeff Re γ} decay against J·J growth e^{2|a||Im(cγ)|}
function _cycle_tail_dir(μeff, aabs, c)
    rate = θ -> μeff * cos(θ) - 2 * aabs * abs(imag(c * cis(θ)))
    θs = range(-1.55, 1.55, length = 63)
    return θs[argmax(map(rate, θs))]
end

@inline function _cycle_rescale(w, logscale)
    aw = abs(w)
    return iszero(aw) ? zero(ComplexF64) : (w / aw) * exp(log(aw) + logscale)
end

function _cycle_harmonic(n, scaledUcov, μeff, ω, Ω, kz, a;
        rtol = 1.0e-7, maxevals = 10^5)
    aabs = abs(a)
    iszero(aabs) && abs(n) > 1 && return zero(AType)   # 𝓣ₙ ≡ 0: keep quadgk off exact zeros
    N = n * Ω
    γ0, D = _cycle_endpoint_D(N, ω, kz)
    c = sqrt(complex(1 - (ω / kz)^2))
    u0 = (γ0 * ω - N) / kz
    e = cis(_cycle_tail_dir(μeff, aabs, c))
    val = QuadGK.quadgk(0.0, Inf; rtol, norm = NORM, maxevals) do t
        δ = t * e
        γ = γ0 + δ
        u = (γ * ω - N) / kz
        q = sqrt(_cycle_q2(γ, γ0, u0, ω, kz; D, δ))   # 𝓣ₙ even in qₙ: sqrt branch irrelevant
        z = a * q
        pref = scaledUcov(γ, u, 2 * abs(imag(z)))
        (pref * e) .* _T_n_bare_x(n, z, u, q)
    end[1]
    return (-2π * im / abs(kz)) .* val
end

# decay-scale estimate μeff ≈ −∂_γ log f from the ∂_γf ratio on the energy axis
function _mueff_estimate(denergy)
    r = abs(first(denergy(1.5, 0.0))) / abs(first(denergy(2.5, 0.0)))
    return clamp(log(r), 0.05, 1.0e3)
end

# Same units as _coupled_contribution: (Pi2/ω²)·X, NOT divided by density.
function _cycle_shell_exponent(μeff, ω, Ω, kz, kperp)
    δ = abs(real(ω)) - abs(kz)
    δ ≤ 0 && return -Inf
    d = abs(complex(δ, imag(ω)))
    iszero(d) && return Inf
    c = abs(sqrt(complex(1 - (ω / kz)^2)))
    x = abs(kperp) * c / d
    η = iszero(x) ? Inf : x < 1 ? log((1 + sqrt(1 - x^2)) / x) - sqrt(1 - x^2) : 0.0
    return μeff * abs(Ω) * δ / d^2 - 2η
end

function _cycle_contribution(c::PreparedVDF, s, ω, k; mueff = nothing, scaledUcov = nothing,
        quad = BoxQuad(_GL24, _GL32), rtol = 1.0e-6, kwargs...)
    box = c.vdf
    denergy = box.denergy
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    a = kperp / Ω
    if isnothing(scaledUcov)
        scaledUcov = (γ, u, σ) -> begin
            dg = denergy(γ, u)
            _cycle_rescale(ω * dg[1] + kz * dg[2], σ)
        end
    end
    μeff = @something mueff _mueff_estimate(denergy)
    nmax = nmax_bessel(a^2 * box.perp[2]^2 / 2)
    f = n -> _harmonic_rel(n, box, ω, Ω, kz, a, quad; landau = false) .+
        _cycle_harmonic(n, scaledUcov, μeff, ω, Ω, kz, a; rtol = 0.1 * rtol)
    # Optimal truncation, NOT plain convergence: shell n pits monodromy
    # e^{μeff·n|Ω|(Reω−k∥)/|k∥−ω|²} against Bessel decay e^{−2n·η(k⊥|c|/|k∥−ω|)}.
    # Term-wise continuation converges iff 2η > μeff|Ω|(Reω−k∥)/|k∥−ω|².
    X_T = f(0)
    divergent = !(_cycle_shell_exponent(μeff, ω, Ω, kz, kperp) < 0)
    # Low orders can rise before Bessel turnover; only sustained post-burn-in growth identifies resurgence.
    burnin = min(nmax, max(3, ceil(Int, abs(a) * box.perp[2]) + 2))
    sprev = Inf
    sbest = Inf
    Xbest = X_T
    nbest = 0
    ngrowing = 0
    for n in 1:nmax
        shell = f(n) + f(-n)
        sn = NORM(shell)
        if !isfinite(sn)
            nbest > 0 || throw(DomainError(sn, "residue-cycle harmonic shell is non-finite before a truncation minimum"))
            X_T = Xbest
            break
        end
        X_T += shell
        if sn < sbest
            sbest, Xbest, nbest = sn, X_T, n
            ngrowing = 0
        elseif divergent && n ≥ burnin
            ngrowing = sn > sprev ? ngrowing + 1 : 0
        end
        if divergent && n ≥ burnin && ngrowing ≥ 3
            X_T = Xbest
            floor_rel = sbest / NORM(X_T)
            floor_rel > 100 * rtol && @warn "residue-cycle harmonic sum is asymptotic here (divergent monodromy tail): optimally truncated at |n| = $nbest with ambiguity floor ≈ $(round(floor_rel, sigdigits = 2)) — outside the convergence domain the continued sheet is defined only to this accuracy" maxlog = 1
            break
        end
        !divergent && n ≥ burnin && sn ≤ rtol * NORM(X_T) && break
        sprev = sn
    end
    X = _antisymmat(2π * X_T) .+ _ee33(c.cache.bernstein33)
    return SMatrix{3, 3, ComplexF64}((s.Pi2 / ω^2) * X)
end
