"""
    Muller(; atol=1e-10, maxiter=100)

Complex-valued Muller's method for derivative-free root polishing.
"""
Base.@kwdef struct Muller
    atol::Float64 = 1.0e-10
    maxiter::Int = 100
end

seedsolver(prob, alg::Muller; kw...) = _seedsolver(prob, alg, _options(kw))
findroot(alg::Muller, D, ω0) = muller(erase_cf(D), ω0; alg.atol, alg.maxiter)

function polish!(f, ωs, alg = Muller())
    n = 0
    for i in eachindex(ωs)
        ωs[i], ni, _ = muller(f, ωs[i]; alg.atol, alg.maxiter)
        n += ni
    end
    return ωs, n
end
