include("Gaussian.jl")
include("ColdVDF.jl")
include("separable.jl")
include("Maxwellian.jl")
include("RingBeam.jl")
include("MaxwellJuttner.jl")
include("GridVDF.jl")
include("SeparableVDF.jl")
include("CoupledVDF.jl")

# Trait wiring
Regime(::AbstractVDF) = NonRelativistic()

