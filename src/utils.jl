const RANGE_DOC = """
`para`/`perp` are the `(lo, hi)` velocity-space integration bounds. 
Scalar values are interpreted as: `para=x` ⇒ `(-x, x)`, `perp=x` ⇒ `(0, x)`. 
For `f₀` not centered on zero, provide the explicit tuple."""

@inline _para_range(x::Union{Tuple,AbstractVector}) = (x[1], x[2])
@inline _para_range(x) = (-x, x)
@inline _perp_range(x::Union{Tuple,AbstractVector}) = (x[1], x[2])
@inline _perp_range(x) = (zero(x), x)