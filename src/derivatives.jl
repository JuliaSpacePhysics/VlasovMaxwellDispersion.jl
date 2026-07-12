# ForwardDiff rejects a complex Dual perturbation;
# a REAL Dual flowing through holomorphic complex arithmetic still yields f'(z)
using ForwardDiff: Dual, value, partials, Partials

struct HoloTag end

@inline _holo_t(t, ::Real) = t
@inline _holo_t(t, x) = complex(t, imag(x))

# value/k-th partial, uniform over Dual and Complex{Dual} results
@inline _val(d::Dual) = value(d)
@inline _val(d::Complex) = complex(value(real(d)), value(imag(d)))
@inline _part(d::Dual, k) = partials(d, k)
@inline _part(d::Complex, k) = complex(partials(real(d), k), partials(imag(d), k))

@inline _dwrt(f, x) = last(_val_dwrt(f, x))
# (f(x), f'(x)) from ONE Dual pass
@inline function _val_dwrt(f, x)
    d = f(_holo_t(Dual{HoloTag}(real(x), one(real(x))), x))
    return _val(d), _part(d, 1)
end

# (∂₁f, ∂₂f) of a 2-arg holomorphic f from ONE pass
@inline function _grad2(f, x, y)
    sx, sy = one(real(x)), one(real(y))
    d = f(
        _holo_t(Dual{HoloTag}(real(x), (sx, zero(sx))), x),
        _holo_t(Dual{HoloTag}(real(y), (zero(sy), sy)), y)
    )
    return _part(d, 1), _part(d, 2)
end

# erf family and gamma functions ship only concrete-Complex methods → MethodError at Complex{Dual}
# Rule: evaluate f at the Complex value, chain the holomorphic derivative: ∂f(z)/∂sₖ = f′(z)·∂z/∂sₖ.
using SpecialFunctions: digamma

@inline function _holo_chain(f, f′, z::Complex{<:Dual{T, V, N}}) where {T, V, N}
    z0 = _val(z)
    w0, d = f(z0), f′(z0)
    wp = ntuple(k -> d * _part(z, k), Val(N))
    return complex(
        Dual{T}(real(w0), Partials(map(real, wp))),
        Dual{T}(imag(w0), Partials(map(imag, wp)))
    )
end

# HoloTag-only signature avoids piracy
SpecialFunctions.erf(z::Complex{<:Dual{HoloTag}}) = _holo_chain(erf, w -> 2 / sqrt(π) * exp(-w^2), z)
SpecialFunctions.erfc(z::Complex{<:Dual{HoloTag}}) = _holo_chain(erfc, w -> -2 / sqrt(π) * exp(-w^2), z)
SpecialFunctions.erfcx(z::Complex{<:Dual{HoloTag}}) = _holo_chain(erfcx, w -> 2w * erfcx(w) - 2 / sqrt(π), z)
SpecialFunctions.erfi(z::Complex{<:Dual{HoloTag}}) = _holo_chain(erfi, w -> 2 / sqrt(π) * exp(w^2), z)
SpecialFunctions.dawson(z::Complex{<:Dual{HoloTag}}) = _holo_chain(dawson, w -> 1 - 2w * dawson(w), z)
SpecialFunctions.gamma(z::Complex{<:Dual{HoloTag}}) = _holo_chain(SF.gamma, w -> SF.gamma(w) * SF.digamma(w), z)
Gamma.gamma(z::Complex{<:Dual{HoloTag}}) = _holo_chain(Gamma.gamma, w -> Gamma.gamma(w) * Gamma.digamma(w), z)
