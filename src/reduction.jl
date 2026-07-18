"""
    TensorReduction

Tensor → scalar reduction whose zeros a dispersion solve chases.

- `:det` — `det 𝒟` - default.
- `:L` / `:R` — circular factors at parallel `k` (`k⊥ = 0`), where
  `det = L·R·P`. `:L` resonates with positive charges (`ω → +Ω`); `:R` mirror.
- `:P` — `k̂ᵀ𝒟k̂`: at `k⊥ = 0` the exact parallel electrostatic factor;
  obliquely the electrostatic approximation.
- `:O` / `:X` — perpendicular (`k∥ = 0`) factors `det = O·X`; exact only when
  every species' `f₀` is even in `p∥` ([`parallel_even`](@ref)).

See docs page *Mode reduction* for physics and derivations.
"""
abstract type TensorReduction end
struct FullDet <: TensorReduction end
struct Circular <: TensorReduction
    σ::Int
    Circular(σ) = abs(σ) == 1 ? new(σ) :
                  throw(ArgumentError("Circular polarization σ must be ±1, got $σ"))
end
struct Longitudinal <: TensorReduction end
struct Ordinary <: TensorReduction end
struct Extraordinary <: TensorReduction end

TensorReduction(m::TensorReduction) = m
@inline TensorReduction(m::Symbol) =
    m === :det ? FullDet() :
    m === :L ? Circular(+1) :
    m === :R ? Circular(-1) :
    m === :P ? Longitudinal() :
    m === :O ? Ordinary() :
    m === :X ? Extraordinary() :
    throw(ArgumentError("unknown mode $m; use :det, :L, :R, :P, :O, :X or a TensorReduction"))

# Validity domain: a geometry condition (checked on the same parameter grid
# the solvers will sample) and, for the perpendicular factors, a plasma
# symmetry condition. The k-component that must vanish is compared against the
# other so a rounded axis point (e.g. k⊥ = k sin(θ), θ ≈ π) still qualifies.
_negligible(x, ref) = abs(x) ≤ 4 * eps(one(x)) * abs(ref)
_valid_at(::TensorReduction, k) = true
_valid_at(::Circular, k) = _negligible(perp(k), para(k))
_valid_at(::Union{Ordinary,Extraordinary}, k) = _negligible(para(k), perp(k))

_valid_plasma(::TensorReduction, plasma) = true
_valid_plasma(::Union{Ordinary,Extraordinary}, plasma) =
    all(s -> parallel_even(s.vdf), NormalizedPlasma(plasma).species)

_domain(::TensorReduction) = "a restricted k domain"
_domain(::Circular) = "exactly parallel k (k⊥ = 0)"
_domain(::Union{Ordinary,Extraordinary}) = "exactly perpendicular k (k∥ = 0)"

function check_reduction(mode, plasma, geometry)
    m = TensorReduction(mode)
    m isa FullDet && return m
    _valid_at(m, geometry) || throw(
        ArgumentError("$(typeof(m)) is defined only at $(_domain(m)), which the geometry leaves")
    )
    _valid_plasma(m, plasma) || throw(
        ArgumentError(
            "$(typeof(m)) needs every f₀ even in p∥, which a species cannot certify " *
            "(drift, or a data-driven VDF); declare " *
            "`VlasovMaxwellDispersion.parallel_even(::YourVDF) = true` to opt in, or use mode = :det"
        )
    )
    return m
end

# Reduction of the tensor to the tracked scalar.
# Circular convention: with e^{-iωt} and B₀ ∥ ẑ, 
# the σ=+1 (L) factor — resonant with positive charges
@inline (m::TensorReduction)(M, k) = m(M)
@inline (::FullDet)(M) = det(M)
@inline (c::Circular)(M) = M[1, 1] + c.σ * im * M[1, 2]
@inline function (::Longitudinal)(M, k)
    kv = vec3(k)
    n2 = sum(abs2, kv)
    return iszero(n2) ? M[3, 3] : dot(kv, M, kv) / n2
end
# Perpendicular (k∥=0, f₀ even in p∥) block factors: det = O·X there.
@inline (::Ordinary)(M) = M[3, 3]
@inline (::Extraordinary)(M) = M[1, 1] * M[2, 2] - M[1, 2] * M[2, 1]

# A factor scales like a product of tensor eigenvalues — 1 of 3 for the 1×1
# factors, 2 of 3 for the Extraordinary 2×2 block; shrinks the det-based
# residual scale accordingly.
_scalepow(::TensorReduction) = 1 / 3
_scalepow(::Extraordinary) = 2 / 3
_scalepow(::FullDet) = 1.0