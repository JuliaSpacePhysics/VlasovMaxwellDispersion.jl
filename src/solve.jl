function residual(prob::DispersionProblem, ω, k = prob.k)
    p, u = normalized(remake(prob; k))
    opts, _ = _split_options(prob.kwargs)
    return residual(p.plasma, ω / u.omega, p.k; closure = _closure(opts.backend))
end

include("solver/muller.jl")
include("solver/GRPF.jl")
include("solver/AAA.jl")
include("solver/survey.jl")

