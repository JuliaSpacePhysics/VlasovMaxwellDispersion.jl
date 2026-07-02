# Holomorphic-function autodiff bridge.
#
# derivatives.jl differentiates a VDF by seeding a *real* ForwardDiff `Dual` and letting it ride
# the Landau contour as the real part of a `Complex{Dual}` (see `_holo_t`). Elementary ops
# (+ - * / ^ exp cos …) already compose through `Complex{Dual}`, so an analytic VDF built from
# them differentiates for free. But the Faddeeva/error-function family (`erf`, `erfc`, `erfcx`,
# …) and `gamma` ship only `Complex{Float64}` methods — evaluating them at `Complex{Dual}` is a
# `MethodError`, which is why an `erfc`-based VDF (e.g. the electron deficit) needed a
# hand-supplied derivative.
#
# Each rule below evaluates `f` at the `Complex{Float64}` value and chain-rules its *holomorphic*
# derivative `f′` through the dual partials: for holomorphic f, ∂/∂sₖ f(z) = f′(z)·∂z/∂sₖ. Cost is
# one `f` + one `f′` call — no truncation, no basis expansion — so an arbitrary analytic VDF made
# of these primitives now differentiates automatically at the same speed as a hand-coded gradient.
#
# Dispatch is restricted to VMD's own `Dual{HoloTag}` (derivatives.jl): the tag is a VMD type in
# the signature, so these are not type-piracy, and no other package's `Complex{Dual}` is affected.
using ForwardDiff: Dual, value, partials, Partials
using SpecialFunctions: digamma

@inline function _holo_chain(f, f′, z::Complex{<:Dual{HoloTag, V, N}}) where {V, N}
    zr, zi = reim(z)
    z0 = complex(value(zr), value(zi))
    w0 = f(z0)
    d  = f′(z0)                                    # holomorphic derivative at z0 (Complex)
    wp = ntuple(k -> d * complex(partials(zr, k), partials(zi, k)), Val(N))
    return complex(Dual{HoloTag}(real(w0), Partials(map(real, wp))),
                   Dual{HoloTag}(imag(w0), Partials(map(imag, wp))))
end

# f′ expressed via the same primitives (all Complex-valued here). √π as a Float is fine — z is Complex.
SpecialFunctions.erf(z::Complex{<:Dual{HoloTag}})    = _holo_chain(erf,    w -> 2 / sqrt(π) * exp(-w^2), z)
SpecialFunctions.erfc(z::Complex{<:Dual{HoloTag}})   = _holo_chain(erfc,   w -> -2 / sqrt(π) * exp(-w^2), z)
SpecialFunctions.erfcx(z::Complex{<:Dual{HoloTag}})  = _holo_chain(erfcx,  w -> 2w * erfcx(w) - 2 / sqrt(π), z)
SpecialFunctions.erfi(z::Complex{<:Dual{HoloTag}})   = _holo_chain(erfi,   w -> 2 / sqrt(π) * exp(w^2), z)
SpecialFunctions.dawson(z::Complex{<:Dual{HoloTag}}) = _holo_chain(dawson, w -> 1 - 2w * dawson(w), z)
SpecialFunctions.gamma(z::Complex{<:Dual{HoloTag}})  = _holo_chain(gamma,  w -> gamma(w) * digamma(w), z)

# Construction-time guard for autodiff-only arbitrary VDFs. `probe` runs the autodiff path once
# at a complex (Landau-contour) point; a `MethodError` means some special function in `f` has no
# `Complex{Dual}` method, so we raise actionable guidance instead of a cryptic failure at solve
# time. Any other probe failure (e.g. an overflow at the probe point) is ignored so it can't block
# construction — the real error, if any, still surfaces during `solve`.
function _assert_holo_diff(probe)
    try
        probe()
    catch err
        err isa MethodError || return nothing
        throw(ArgumentError(
            "VDF is not auto-differentiable at a complex (Landau-contour) argument: " *
            "`$(err.f)` has no `Complex{Dual}` method. Either pass an explicit gradient " *
            "(`dfpara`/`dfperp`/`dgrad`), or build the VDF only from primitives the holomorphic-AD " *
            "bridge covers — elementary ops (+ - * / ^ exp log trig) and erf/erfc/erfcx/erfi/dawson/gamma."))
    end
    return nothing
end
