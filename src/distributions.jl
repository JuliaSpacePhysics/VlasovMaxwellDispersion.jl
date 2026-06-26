include("distributions/ColdVDF.jl")
include("distributions/Maxwellian.jl")
include("distributions/RingBeam.jl")
include("distributions/MaxwellJuttner.jl")
include("distributions/GridVDF.jl")
include("distributions/SeparableVDF.jl")
include("distributions/CoupledVDF.jl")

# Trait wiring
Regime(::AbstractVDF) = NonRelativistic()

