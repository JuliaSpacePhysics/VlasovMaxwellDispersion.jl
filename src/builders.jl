# Physical `Plasma` → this solver's dimensionless form

"""
    scales(plasma::Plasma; ref = plasma.populations[1].particle)

Dimensionless scales of a physical plasma: speeds in `c`, frequencies in `Ω_ref`, lengths in
`c/Ω_ref`. `ref` is a `Particle`, whose gyrofrequency is `Ω_ref`, or a frequency [rad/s].

- `vA`, `B0` [T], `Omega_ref` [rad/s]
- `populations`: per population `(; wc, wp, d, vth_perp, vth_para, rho)` with signed
  `wc = Ω_s/Ω_ref`, magnitude `wp = ω_ps/Ω_ref`, inertial length `d = 1/wp`, and, for PlasmaPopulations
  models, thermal speeds and gyroradius `rho = vth_perp/|wc|` (else `nothing`)
"""
function scales(plasma::Plasma; ref = plasma.populations[1].particle)
    B0 = norm(plasma.B0)
    Ωref = ref isa Particle ? abs(gyrofrequency(ref, B0)) : si_frequency(ref)
    sp = map(plasma.populations) do s
        wc = gyrofrequency(s.particle, B0) / Ωref
        wp = plasma_frequency(s.particle, s.n) / Ωref
        vth_perp, vth_para = _vth_perp_para(s.vdf)
        rho = isnothing(vth_perp) ? nothing : vth_perp / abs(wc)
        (; wc, wp, d = 1 / wp, vth_perp, vth_para, rho)
    end
    return (; vA = alfven_speed(plasma) / C_SI, B0, Omega_ref = Ωref, populations = sp)
end

_vth_perp_para(d::AbstractVDF) = vth_perp_para(d) ./ C_SI
_vth_perp_para(_) = (nothing, nothing)

function NormalizedPlasma(plasma::Plasma; ref = plasma.populations[1].particle)
    sc = scales(plasma; ref)
    sps = map(plasma.populations, sc.populations) do s, ss
        _normalized_population(ss.wc, sign(s.n) * ss.wp^2, s.vdf, C_SI)
    end
    return NormalizedPlasma(sps, sc)
end

prepare(p::Plasma, args...; kw...) = prepare(NormalizedPlasma(p), args...; kw...)
