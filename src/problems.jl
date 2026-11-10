# This solver's `solve` keywords, stored on a `DispersionProblem` or given to `solve`.
const _OPTION_DEFAULTS = (; mode = FullDet(), backend = FixedNodeEval())

# `kw` split into resolved options and the keywords an algorithm passes on to a wrapped solver.
function _split_options(kw)
    kw = NamedTuple(kw)
    rest = Base.structdiff(kw, _OPTION_DEFAULTS)
    o = merge(_OPTION_DEFAULTS, Base.structdiff(kw, rest))
    return (; mode = TensorReduction(o.mode), o.backend), rest
end

# Options of an algorithm that wraps no other solver, so a misspelled or foreign keyword is an error.
function _options(kw)
    opts, rest = _split_options(kw)
    isempty(rest) ||
        throw(ArgumentError("unknown keywords $(keys(rest)); this solver takes $(keys(_OPTION_DEFAULTS))"))
    return opts
end

# The normalized problem with its populations prepared for `opts`. Checking `mode` at every
# wavevector here fails a solve before its first root rather than mid-path or inside a threaded survey.
function _prepare(prob, opts)
    p = remake(prob; plasma = prepare(prob.plasma; opts.backend))
    foreach(k -> check_reduction(opts.mode, p.plasma, k), wavevectors(p.k))
    return p
end
