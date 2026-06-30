# ForwardDiff rejects a complex Dual perturbation;
# a REAL Dual flowing through holomorphic complex arithmetic still yields f'(z)
using ForwardDiff: derivative, Dual, value, partials

# VMD-owned ForwardDiff tag. Every differentiation seed below carries it, so the
# holomorphic special-function rules (holomorphic_ad.jl) dispatch on `Dual{HoloTag}` —
# a VMD type in the signature, hence not type piracy — and never touch anyone else's Duals.
struct HoloTag end

@inline _holo_t(t, ::Real) = t # performance optimization so perpendicular (real-argument) moments stay real-typed
@inline _holo_t(t, x) = complex(t, imag(x))
@inline _seed1(x) = Dual{HoloTag}(real(x), one(real(x)))
@inline _dwrt(f, x) = last(_splitdual(f(_holo_t(_seed1(x), x))))

# (f(x), f'(x)) from ONE Dual pass.
@inline _val_dwrt(f, x) = _splitdual(f(_holo_t(_seed1(x), x)))
@inline _splitdual(d::Dual) = (value(d), partials(d, 1))
@inline function _splitdual(d::Complex)
    dr, di = reim(d)
    return complex(value(dr), value(di)), complex(partials(dr, 1), partials(di, 1))
end

# (∂₁f, ∂₂f) of a 2-arg holomorphic f from ONE passs
@inline function _grad2(f, x, y)
    sx, sy = one(real(x)), one(real(y))
    dx = _holo_t(Dual{HoloTag}(real(x), (sx, zero(sx))), x)
    dy = _holo_t(Dual{HoloTag}(real(y), (zero(sy), sy)), y)
    return _split2(f(dx, dy))
end
@inline _split2(d::Dual) = (partials(d, 1), partials(d, 2))
@inline function _split2(d::Complex)
    dr, di = reim(d)
    return complex(partials(dr, 1), partials(di, 1)), complex(partials(dr, 2), partials(di, 2))
end
