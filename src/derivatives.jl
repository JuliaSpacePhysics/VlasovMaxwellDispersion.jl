# ForwardDiff rejects a complex Dual perturbation; 
# a REAL Dual flowing through holomorphic complex arithmetic still yields f'(z)
using ForwardDiff: derivative

@inline _holo_t(t, ::Real) = t # performance optimization so perpendicular (real-argument) moments stay real-typed
@inline _holo_t(t, x) = complex(t, imag(x))
@inline _dwrt1(f, u, v) = derivative(t -> f(_holo_t(t, u), v), real(u))  # ∂/∂(1st)
@inline _dwrt2(f, u, v) = derivative(t -> f(u, _holo_t(t, v)), real(v))  # ∂/∂(2nd)
@inline _dwrt(f, x) = derivative(t -> f(_holo_t(t, x)), real(x))         # 1-arg
