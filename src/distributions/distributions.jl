include("Gaussian.jl")
include("cold.jl")
include("separable.jl")
include("Maxwellian.jl")
include("RingBeam.jl")
include("MaxwellJuttner.jl")
include("SplineVDF.jl")
include("ReducedVDF.jl")
include("relativistic.jl")
include("CoupledVDF.jl")
include("residue_cycles.jl")
include("LowRankVDF.jl")
include("Kappa.jl")
include("BiKappa.jl")

"""
    solver_vdf(vdf, vc)

`vdf` in this solver's representation. PlasmaPopulations models are rescaled by `vc`, the speed
of light in their speed unit; a [`SolverVDF`](@ref) passes through.
"""
solver_vdf(d::SolverVDF, vc) = d
solver_vdf(d, vc) = throw(ArgumentError("$(nameof(typeof(d))) has no representation in this solver"))

_inc(v, vc) = v / vc

_check_nodrift(d) = iszero(d.vd_para) ||
    throw(ArgumentError("$(nameof(typeof(d))) with a parallel drift `vd_para` is not supported"))

function solver_vdf(d::BiMaxwellian, vc)
    vperp, vpara, vd = map(v -> _inc(v, vc), (d.vth_perp, d.vth_para, d.vd_para))
    return Gaussian(vperp) ⊗ Gaussian(vpara, vd)
end

function solver_vdf(d::GaussianRing, vc)
    w_perp, vpara, vd, vr = map(v -> _inc(v, vc), (d.w_perp, d.vth_para, d.vd_para, d.vr_perp))
    perp = iszero(vr) ? Gaussian(w_perp) : Gaussian(w_perp, vr)
    return perp ⊗ Gaussian(vpara, vd)
end

function solver_vdf(d::BiKappa, vc)
    _check_nodrift(d)
    κ = d.kappa
    a_perp, a_para = (κ - 3 / 2) .* _inc.((d.vth_perp, d.vth_para), vc) .^ 2
    return CoupledKappa(isinteger(κ) ? Int(κ) : κ, promote(a_para, a_perp)...)
end

function solver_vdf(d::ProductBiKappa, vc)
    _check_nodrift(d)
    iszero(d.loss_cone) || throw(ArgumentError("ProductBiKappa with a loss cone is not supported"))
    θperp, θpara = _inc.(kappa_widths(d), vc)
    return Kappa(θperp, d.kappa_perp) ⊗ Kappa(θpara, d.kappa_para)
end

solver_vdf(d::Cold, vc) = (_check_nodrift(d); d)
solver_vdf(d::MaxwellJuttner, vc) = Juttner(inv(float(d.theta)))
solver_vdf(d::GridVDF, vc) = SplineVDF(d.v_perp ./ vc, d.v_para ./ vc, d.f)
