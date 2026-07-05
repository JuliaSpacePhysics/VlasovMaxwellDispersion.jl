# Prototype: infinite-bounds non-rel CoupledVDF harmonic path.
# Replaces finite-box Plemelj subtraction (g-gζ, log((U-ζ)/(L-ζ))) — divergent on ℝ —
# with Lorentzian-kernel subtraction:
#   ∫ g/(u-ζ) du = ∫ [g - gζ·h(u)]/(u-ζ) du + gζ·C(ζ),  h = W²/(W²+(u-Re ζ)²)
#   C(ζ) = iπW/(W+Im ζ)          Im ζ ≥ 0
#        = 2πi − iπW/(W−Im ζ)    Im ζ < 0 (Landau-crossed: real ζ always "inside" ℝ)
using VlasovMaxwellDispersion
using VlasovMaxwellDispersion: _In_block, _perp_Bessel_bilinears!, _relsize, _subtract_safe, _antisymmat,
    nmax_bessel, AType, para, perp, contribution, NORM
using QuadGK, StaticArrays, Chairmarks
using Bumper: @no_escape, @alloc

_C_inf(ζ, W) = imag(ζ) >= 0 ? im * π * W / (W + imag(ζ)) : 2π * im - im * π * W / (W - imag(ζ))

function _coupled_perp_inf(v, ns, ζs, d, ω, Ω, kz, a, W; kw...)
    g5(u) = begin
        q, p = d.dgrad(v, u)
        SVector(q, u * q, u^2 * q, p, u * p)
    end
    invkz = -1 / kz
    nb = length(ns)
    gscale = maximum(ζ -> _relsize(g5(real(ζ))), ζs)
    W2 = W^2
    return @no_escape begin
        gζs = @alloc(SVector{5, eltype(ζs)}, nb)
        near = @alloc(Bool, nb)
        @inbounds for i in 1:nb
            gζs[i] = g5(ζs[i])
            near[i] = _subtract_safe(gζs[i], gscale)
        end
        b2s = @alloc(SVector{6, typeof(a * v)}, nb)
        _perp_Bessel_bilinears!(b2s, a, v)
        reg = QuadGK.quadgk(-Inf, Inf; kw...) do u
            g = g5(u)
            acc = zero(AType)
            @inbounds for i in eachindex(ns)
                c = invkz / (u - ζs[i])
                h = W2 / (W2 + (u - real(ζs[i]))^2)
                acc += _In_block(near[i] ? g - h * gζs[i] : g, c, b2s[i], v, ω, kz, ns[i] * Ω)
            end
            acc
        end[1]
        logacc = zero(AType)
        @inbounds for i in eachindex(ns)
            pc = near[i] ? gζs[i] .* _C_inf(ζs[i], W) :
                (imag(ζs[i]) < 0 ? gζs[i] .* (2π * im) : zero(gζs[i]))
            logacc += _In_block(pc, invkz, b2s[i], v, ω, kz, ns[i] * Ω)
        end
        reg + logacc
    end
end

# Gaussian kernel: h = exp(-((u-Re ζ)/W)²), C = √π·Z-type closed form via erfcx:
#   Im ζ ≥ 0: iπ·erfcx(Im ζ/W);  Im ζ < 0: 2πi − iπ·erfcx(−Im ζ/W)
using SpecialFunctions: erfcx
_C_gauss(ζ, W) = imag(ζ) >= 0 ? im * π * erfcx(imag(ζ) / W) : 2π * im - im * π * erfcx(-imag(ζ) / W)

function _coupled_perp_inf2(v, ns, ζs, d, ω, Ω, kz, a, W; kw...)
    g5(u) = begin
        q, p = d.dgrad(v, u)
        SVector(q, u * q, u^2 * q, p, u * p)
    end
    invkz = -1 / kz
    nb = length(ns)
    gscale = maximum(ζ -> _relsize(g5(real(ζ))), ζs)
    return @no_escape begin
        gζs = @alloc(SVector{5, eltype(ζs)}, nb)
        near = @alloc(Bool, nb)
        @inbounds for i in 1:nb
            gζs[i] = g5(ζs[i])
            near[i] = _subtract_safe(gζs[i], gscale)
        end
        b2s = @alloc(SVector{6, typeof(a * v)}, nb)
        _perp_Bessel_bilinears!(b2s, a, v)
        reg = QuadGK.quadgk(-Inf, Inf; kw...) do u
            g = g5(u)
            acc = zero(AType)
            @inbounds for i in eachindex(ns)
                c = invkz / (u - ζs[i])
                h = exp(-((u - real(ζs[i])) / W)^2)
                acc += _In_block(near[i] ? g - h * gζs[i] : g, c, b2s[i], v, ω, kz, ns[i] * Ω)
            end
            acc
        end[1]
        logacc = zero(AType)
        @inbounds for i in eachindex(ns)
            pc = near[i] ? gζs[i] .* _C_gauss(ζs[i], W) :
                (imag(ζs[i]) < 0 ? gζs[i] .* (2π * im) : zero(gζs[i]))
            logacc += _In_block(pc, invkz, b2s[i], v, ω, kz, ns[i] * Ω)
        end
        reg + logacc
    end
end

function _contribution_infcore(perpfun::F, d, s, ω, k, W, rtol, norm) where {F}
    ω = complex(float(ω))
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    a = kperp / Ω
    p⊥²_mean = 2π * QuadGK.quadgk(
        v -> v^3 * QuadGK.quadgk(u -> d.f0(v, u), -Inf, Inf; rtol = 1.0e-3)[1],
        0.0, Inf; rtol = 1.0e-3
    )[1] / d.n
    nmax = nmax_bessel(a^2 * abs(p⊥²_mean) / 2)
    ns = (-nmax):nmax
    ζs = [(ω - n * Ω) / kz for n in ns]
    X = QuadGK.quadgk(0.0, Inf; rtol, norm) do v
        perpfun(v, ns, ζs, d, ω, Ω, kz, a, W; norm, rtol)
    end[1]
    return (s.Pi2 / ω^2) * _antisymmat(X) / d.n
end

contribution_inf(d, s, ω, k; W = 1.0, rtol = 1.0e-6, norm = NORM) =
    _contribution_infcore(_coupled_perp_inf, d, s, ω, k, W, rtol, norm)
contribution_inf2(d, s, ω, k; W = 1.0, rtol = 1.0e-6, norm = NORM) =
    _contribution_infcore(_coupled_perp_inf2, d, s, ω, k, W, rtol, norm)
"proto loaded"

# WINNER — manual map ℝ→(−1,1), u = S·t/(1−t²), then CONSTANT Plemelj subtraction in t.
# Residue of jac·g/(φ(t)−ζ) at t* = φ⁻¹(ζ) is exactly g(ζ) (Jacobian cancels), so the
# box machinery ports verbatim: subtract gζ/(t−t*), add gζ·[log((1−t*)/(−1−t*)) + 2πi if
# Im ζ<0]. Ghost preimage t̃ = −1/t* lies outside (−1,1) and only matters where g≈0.
# S ≈ 4·√⟨p∥²⟩ (broad optimum 2–8×); outer perp stays builtin quadgk(0,∞) — scaled/mapped
# outer variants regressed (ring 243 vs 164 ms; Maxwellian 11.8 vs 9.2 ms).
_tstar(ζS) = 2ζS / (1 + sqrt(1 + 4ζS^2))

function _coupled_perp_map(v, ns, ζs, d, ω, Ω, kz, a, S; kw...)
    g5(u) = begin
        q, p = d.dgrad(v, u)
        SVector(q, u * q, u^2 * q, p, u * p)
    end
    invkz = -1 / kz
    nb = length(ns)
    gscale = maximum(ζ -> _relsize(g5(real(ζ))), ζs)
    return @no_escape begin
        gζs = @alloc(SVector{5, eltype(ζs)}, nb)
        near = @alloc(Bool, nb)
        ts = @alloc(eltype(ζs), nb)
        @inbounds for i in 1:nb
            gζs[i] = g5(ζs[i])
            near[i] = _subtract_safe(gζs[i], gscale)
            ts[i] = _tstar(ζs[i] / S)
        end
        b2s = @alloc(SVector{6, typeof(a * v)}, nb)
        _perp_Bessel_bilinears!(b2s, a, v)
        reg = QuadGK.quadgk(-1.0, 1.0; kw...) do t
            u = S * t / (1 - t^2)
            jac = S * (1 + t^2) / (1 - t^2)^2
            g = g5(u)
            acc = zero(AType)
            @inbounds for i in eachindex(ns)
                cu = jac * invkz / (u - ζs[i])
                m = near[i] ? cu * g - (invkz / (t - ts[i])) * gζs[i] : cu * g
                acc += _In_block(m, 1, b2s[i], v, ω, kz, ns[i] * Ω)
            end
            acc
        end[1]
        logacc = zero(AType)
        @inbounds for i in eachindex(ns)
            crossed = imag(ζs[i]) < 0
            pc = near[i] ? gζs[i] .* (log((1 - ts[i]) / (-1 - ts[i])) + (crossed ? 2π * im : 0)) :
                (crossed ? gζs[i] .* (2π * im) : zero(gζs[i]))
            logacc += _In_block(pc, invkz, b2s[i], v, ω, kz, ns[i] * Ω)
        end
        reg + logacc
    end
end

function contribution_map(d, s, ω, k; cS = 4.0, rtol = 1.0e-6, norm = NORM)
    ω = complex(float(ω))
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    a = kperp / Ω
    pm = 2π * QuadGK.quadgk(v -> v^3 * QuadGK.quadgk(u -> d.f0(v, u), -Inf, Inf; rtol = 1.0e-3)[1], 0.0, Inf; rtol = 1.0e-3)[1] / d.n
    pa = 2π * QuadGK.quadgk(v -> v * QuadGK.quadgk(u -> u^2 * d.f0(v, u), -Inf, Inf; rtol = 1.0e-3)[1], 0.0, Inf; rtol = 1.0e-3)[1] / d.n
    S = cS * sqrt(abs(pa))
    nmax = nmax_bessel(a^2 * abs(pm) / 2)
    ns = (-nmax):nmax
    ζs = [(ω - n * Ω) / kz for n in ns]
    X = QuadGK.quadgk(0.0, Inf; rtol, norm) do v
        _coupled_perp_map(v, ns, ζs, d, ω, Ω, kz, a, S; norm, rtol)
    end[1]
    return (s.Pi2 / ω^2) * _antisymmat(X) / d.n
end

# WINNER v2 — sinh-mapped trapezoid + cotangent pole correction (Kress-type product quadrature).
# Residue calculus on (π/h)cot(πw/h)·G(w)/(ψ(w)−ζ) over the analyticity strip gives, for
# uniform nodes t_j = j·h (alignment REQUIRED: offset grids break the cot phase — bit me),
#   ∫_ℝ g(u)/(u−ζ) du  [Landau-continued]  =  h·Σ_j G(t_j)/(ψ(t_j)−ζ) + π·g(ζ)·(cot(πt*/h) + i)
# with u=ψ(t)=S·sinh(t), G=g(ψ)ψ′, t*=asinh(ζ/S); the residue is g(ζ) exactly (Jacobian
# cancels) and the SAME formula covers Im ω ≷ 0 — the +iπ merges the crossed 2πi into one
# analytic expression, no case split. Error ~e^{−2πd/h} (d = strip width of G): h=0.25→8e-8,
# h=0.2→7e-10, h=0.15→3e-13. Fixed nodes, no adaptivity, all harmonics share the samples.
# S=√⟨p∥²⟩; T=7 (u≈550·S) + f0-amplitude node skip keeps empty tails free.
# Beam/drifted f₀ needs a centered map u = u_c + S·sinh(t), u_c=⟨p∥⟩ (untested).
function _coupled_perp_sinc(v, ns, ζs, d, ω, Ω, kz, a, S; h = 0.2, T = 7.0, ftol = 1.0e-14, kw...)
    g5(u) = begin
        q, p = d.dgrad(v, u)
        SVector(q, u * q, u^2 * q, p, u * p)
    end
    invkz = -1 / kz
    nb = length(ns)
    jmax = floor(Int, T / h)
    return @no_escape begin
        b2s = @alloc(SVector{6, typeof(a * v)}, nb)
        _perp_Bessel_bilinears!(b2s, a, v)
        Is = @alloc(SVector{5, ComplexF64}, nb)
        @inbounds for i in 1:nb
            Is[i] = zero(SVector{5, ComplexF64})
        end
        fmax = abs(d.f0(v, 0.0)) + abs(d.f0(v, S)) + abs(d.f0(v, -S))
        for j in (-jmax):jmax
            t = j * h
            u = S * sinh(t)
            abs(d.f0(v, u)) < ftol * fmax && continue
            w = h * S * cosh(t)
            g = g5(u)
            @inbounds for i in 1:nb
                Is[i] += (w / (u - ζs[i])) * g
            end
        end
        acc = zero(AType)
        @inbounds for i in 1:nb
            tst = asinh(ζs[i] / S)
            gz = g5(ζs[i])
            corr = all(isfinite, gz) ? (π * (cot(π * tst / h) + im)) * gz : zero(gz)
            acc += _In_block(Is[i] + corr, invkz, b2s[i], v, ω, kz, ns[i] * Ω)
        end
        acc
    end
end

function contribution_sinc(d, s, ω, k; h = 0.2, rtol = 1.0e-6, norm = NORM)
    ω = complex(float(ω))
    Ω, kz, kperp = s.Omega, para(k), perp(k)
    a = kperp / Ω
    pm = 2π * QuadGK.quadgk(v -> v^3 * QuadGK.quadgk(u -> d.f0(v, u), -Inf, Inf; rtol = 1.0e-3)[1], 0.0, Inf; rtol = 1.0e-3)[1] / d.n
    pa = 2π * QuadGK.quadgk(v -> v * QuadGK.quadgk(u -> u^2 * d.f0(v, u), -Inf, Inf; rtol = 1.0e-3)[1], 0.0, Inf; rtol = 1.0e-3)[1] / d.n
    S = sqrt(abs(pa))
    nmax = nmax_bessel(a^2 * abs(pm) / 2)
    ns = (-nmax):nmax
    ζs = [(ω - n * Ω) / kz for n in ns]
    X = QuadGK.quadgk(0.0, Inf; rtol, norm) do v
        _coupled_perp_sinc(v, ns, ζs, d, ω, Ω, kz, a, S; h, norm)
    end[1]
    return (s.Pi2 / ω^2) * _antisymmat(X) / d.n
end
