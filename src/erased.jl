# Type-erased user functions to avoid recompilation of entire quadrature
# stacks for every new User-supplied f₀/gradients.
# The solvers actually only use (real node u, complex Landau ζ) as arguments.
# `raw` provides fallback for exotic eltypes (Dual, Float32, BigFloat).
using FunctionWrappers: FunctionWrapper as FW

struct Erased2{R, C}
    rr::FW{R, Tuple{Float64, Float64}}
    rc::FW{C, Tuple{Float64, ComplexF64}}
    raw::Any
end
(f::Erased2)(q::Float64, u::Float64) = f.rr(q, u)
(f::Erased2)(q::Float64, u::ComplexF64) = f.rc(q, u)
(f::Erased2)(q, u) = f.raw(q, u)

struct Erased1{R, C}
    rr::FW{R, Tuple{Float64}}
    rc::FW{C, Tuple{ComplexF64}}
    raw::Any
end
(f::Erased1)(x::Float64) = f.rr(x)
(f::Erased1)(x::ComplexF64) = f.rc(x)
(f::Erased1)(x) = f.raw(x)

# generic shims: always cfunction-compilable even when f has typed methods,
# and normalize 2-value returns to tuples
_shim1(f) = x -> f(x)
_shim2(f) = (q, u) -> f(q, u)
_tup1(f) = x -> ((a, b) = f(x); (a, b))
_tup2(f) = (q, u) -> ((a, b) = f(q, u); (a, b))

const T2 = NTuple{2, Float64}
const CT2 = NTuple{2, ComplexF64}

# erase only on the Float64 lattice (guard = a domain bound); pass through otherwise
erase_f2(f, ::Float64) = f isa Erased2 ? f : Erased2{Float64, ComplexF64}(_shim2(f), _shim2(f), f)
erase_g2(f, ::Float64) = f isa Erased2 ? f : Erased2{T2, CT2}(_tup2(f), _tup2(f), f)
erase_f1(f, ::Float64) = f isa Erased1 ? f : Erased1{Float64, ComplexF64}(_shim1(f), _shim1(f), f)
erase_fd1(f, ::Float64) = f isa Erased1 ? f : Erased1{T2, CT2}(_tup1(f), _tup1(f), f)
erase_f2(f, _) = f
erase_g2(f, _) = f
erase_f1(f, _) = f
erase_fd1(f, _) = f

# # SciML documentation on specialization levels
# https://docs.sciml.ai/DiffEqDocs/stable/features/low_dep/
# https://docs.sciml.ai/SciMLBase/stable/interfaces/Problems/#Specialization-Levels
