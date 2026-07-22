# Relativistic gyrotropic box/peel engine — generic over the VDF
# Math: docs/src/relativistic.typ.

const _GL24 = GaussLegendre(24)
const _GL32 = GaussLegendre(32)

const AType = SVector{6,ComplexF64}

function _warn_damped_superluminal(ω, kz)
    return if imag(ω) < 0 && real(ω)^2 > kz^2
        @warn "damped superluminal ω (|Re ω| > |k∥|): the (p⊥,p∥) integral is not the analytic continuation there (apex branch cut, docs/src/relativistic.typ); use an analytic energy-form VDF with path=:cycles or evaluate at Im ω ≥ 0" maxlog = 1
    end
end

# Covariant momentum numerator 𝒰 = ω∂_γf+k∥∂_uf at (γ,p∥) with w=p⊥, rewritten via
# ∂_γ|_u=(γ/w)∂_⊥, ∂_u|_γ=∂_∥−(u/w)∂_⊥ ⇒ 𝒰 = k∥∂_∥f + (ωγ−k∥u)/w · ∂_⊥f.
@inline function _U_cov(d, u, w, γ, ω, kz)
    dpe, dpa = d.dgrad(w, u)
    return kz * dpa + dpe * (ω * γ - kz * u) / w
end

# Relativistic non-resonant e∥e∥ term without prefactor
# Edge-mapped: γ→q² concentrates nodes near γ=1, p∥→θ half-angle over the resonance ellipse
function _bernstein_rel(d, qs=BoxQuad(_GL24, _GL32))
    γmax = sqrt(1 + max(d.para[1]^2, d.para[2]^2) + d.perp[2]^2)
    acc = quad(qs.outer, 0, 1) do q
        γ = 1 + (γmax - 1) * q^2
        wγ = 2 * (γmax - 1) * q
        umax = sqrt(γ^2 - 1)
        inner = quad(qs.inner, -1, 1) do t
            θ = t * (π / 2)
            u, w = umax .* sincos(θ)
            dpe, dpa = d.dgrad(w, u)
            (π / 2) * ComplexF64(w * u * dpa - u^2 * dpe)
        end
        wγ * inner
    end
    return 2π * acc
end

# Poles farther than this from the real p∥ segment leave the integrand smooth enough
# for plain adaptive quadrature; nearer (or Landau-crossed) poles are peeled.
const _PQ_NEAR = 1.5

# Partial fractions of the rationalized resonance at fixed
# m⊥² = 1+p⊥²: D_n·D̃_n = A·u² + B·u + C, so
#   1/D_n = D̃_n·[c₁/(u−p₁) + c₂/(u−p₂)],  c₁₂ = ∓1/√(B²−4AC).
@inline function _Dn_poles(ω, kz, nΩ, m2)
    A = ω^2 - kz^2
    B = -2 * kz * nΩ
    C = ω^2 * m2 - nΩ^2
    sq = sqrt(B^2 - 4 * A * C)
    abs2(B + sq) < abs2(B - sq) && (sq = -sq)
    return (_home_side((-B - sq) / (2A), ω, kz, m2), -1 / sq),
    (_home_side(2 * C / (-B - sq), ω, kz, m2), 1 / sq)
end

# Exactly-real ω leaves a real pole ON the path: a signed zero nudges it to its home
# side (the Im ω→0⁺ limit, slope dp/dω = γ²/(k∥γ−ωp)) so the boundary-value log in
# `_peel_residue` lands on the correct sheet.
@inline function _home_side(p, ω, kz, m2)
    isfinite(p) || return complex(Inf)
    iszero(imag(p)) || return p
    γp = sqrt(complex(m2 + p^2))
    return complex(real(p), sign(real(γp^2 / (kz * γp - ω * p))) * 0.0)
end

# Residue r = c·W(p) of a peeled pole and its analytic across-box term
# r·[log((hi−p)/(lo−p)) + sgnkz·2πi if Landau-crossed]; (0, 0) when not peeled.
# landau=false suppresses the crossed 2πi terms (near-peel only)
@inline function _peel_residue(p, c, W, γof, ν, lo, hi, sgnkz; landau=true)
    zz = (zero(AType), zero(AType))
    isfinite(p) || return zz
    crossed = landau && ν < 0 && sgnkz * imag(p) < 0 && lo < real(p) < hi
    near = abs(imag(p)) < _PQ_NEAR && lo - _PQ_NEAR < real(p) < hi + _PQ_NEAR
    (crossed || near) || return zz
    γp = γof(p)
    # γ-artifact pole: kz=0,n=0 collapses both rationalized roots onto γ=0 (p=±i√(1+p⊥²)).
    # f₀ is singular (√(1+p⊥²+p∥²)=0, autodiff throws) — short-circuit before evaluating W
    iszero(γp) && return zz
    r = c .* W(p, γp)
    all(isfinite, r) || return zz
    return r, r .* (log((hi - p) / (lo - p)) + (crossed ? sgnkz * 2π * im : 0))
end

# One harmonic of the (p⊥,p∥) box integral: outer Gauss–Legendre in p⊥; per slice
#   ∫ σ𝓣_n/D_n du = ∫ Σᵢ cᵢ·W/(u−pᵢ) du,  W = σ·D̃_n·𝓣_n,  σ = 𝒰·p⊥/γ,
# with peeled poles kept as single fractions (c·W(u)−r)/(u−p) — the split form
# σ𝓣/D − r/(u−p) carries 1/(u−p)² rounding noise near the pole.
function _harmonic_rel(n, d, ω, Ω, kz, a, qs::BoxQuad; landau=true)
    plo, phi = d.para
    qlo, qhi = d.perp
    nΩ = n * Ω
    ν = imag(ω)
    sgnkz = sign(kz)
    return quad(qs.outer, qlo, qhi) do q
        m2 = 1 + q^2
        z = a * q
        γof(u) = sqrt(complex(m2 + u^2))
        σof(u, γ) = _U_cov(d, u, q, γ, ω, kz) * (q / γ)   # slice weight σ = 𝒰 · (p⊥/γ)
        Wof = (u, γ) -> (σof(u, γ) * (ω * γ + kz * u + nΩ)) .* _T_n_bare_x(n, z, u, q)
        (p1, c1), (p2, c2) = _Dn_poles(ω, kz, nΩ, m2)
        r1, lg1 = _peel_residue(p1, c1, Wof, γof, ν, plo, phi, sgnkz; landau)
        r2, lg2 = _peel_residue(p2, c2, Wof, γof, ν, plo, phi, sgnkz; landau)
        reg = quad(qs.inner, plo, phi) do u
            γ = γof(u)
            if isfinite(p1) || isfinite(p2)
                Wu = Wof(u, γ)
                acc = zero(Wu)
                isfinite(p1) && (acc = acc .+ (c1 .* Wu .- r1) ./ (u - p1))
                isfinite(p2) && (acc = acc .+ (c2 .* Wu .- r2) ./ (u - p2))
                acc
            else                       # A=B=0 (ω=±k∥, n=0): quadratic degenerate, no poles
                (σof(u, γ) .* _T_n_bare_x(n, z, u, q)) ./ (ω * γ - kz * u - nΩ)
            end
        end
        reg .+ lg1 .+ lg2
    end
end
