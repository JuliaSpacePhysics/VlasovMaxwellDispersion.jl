# Shared setup for the "transverse heating lowers the EIC current threshold" experiment.
#
# Question: electrostatic ion-cyclotron (EIC) waves are reported in the topside auroral
# ionosphere at field-aligned currents *below* the Kindel & Kennel (1971) current
# threshold. Transversely heated (ring/conic) O⁺ is co-observed with those waves. Does
# that ring free energy lower the current threshold, and by how much?
#
# Model: O⁺ core (fraction 1-δ, isotropic Maxwellian) + O⁺ ring (fraction δ, gyro-ring
# I₀ Maxwellian at ring speed vr) + Maxwellian electrons drifting at u along B₀.
# Ω_ref = Ω_O⁺, speeds in c, k reported as k⊥ρ_i with ρ_i = w_i/Ω_O⁺, w = sqrt(2T/m).
#
# Electrostatic reduction (mode = :P) — the K&K dispersion relation. At these β
# (~1e-4) the EIC/IA band is electrostatic to well below the numbers we quote.

using VlasovMaxwellDispersion
using Printf

const MI_ME = 16 * 1836.15         # O⁺ / electron mass ratio

# --- topside auroral reference point (~1700 km, Freja/FAST altitudes) --------
const B0 = 2.8e-5                  # T
const N0 = 3.0e9                   # m⁻³  (3×10³ cm⁻³)
const TI = 0.3                     # eV   (O⁺ core)
const QE = 1.602176634e-19
const MO = 16 * 1.67262192369e-27
const ME = 9.1093837015e-31
const EPS0 = 8.8541878128e-12
const CLIGHT = 2.99792458e8

const OMEGA_I = QE * B0 / MO                                  # rad/s
const WI = sqrt(2 * TI * QE / MO)                             # O⁺ thermal speed, m/s
const EPS = WI / CLIGHT                                       # w_i/c — the small parameter
const PI2 = N0 * QE^2 / (EPS0 * MO) / OMEGA_I^2               # (ω_pi/Ω_i)²
const RHO_I = WI / OMEGA_I                                    # m

"""Electron thermal speed [m/s] at `tau = Te/Ti`."""
w_e(tau = 1.0) = WI * sqrt(tau * MI_ME)

"""Field-aligned current density [A/m²] from drift `u` in electron thermal speeds."""
current_density(u; tau = 1.0, n = N0) = n * QE * u * w_e(tau)

"""Ring kinetic energy [eV] of the O⁺ ring at ring speed `vr` in core thermal speeds."""
ring_energy(vr) = TI * vr^2        # ½m vr² / e, with w² = 2T/m ⇒ ½m(vr·w)² = vr²·T

# ---------------------------------------------------------------------------
# Plasma builders
# ---------------------------------------------------------------------------

"""O⁺ core + ring + drifting electrons.

`delta` ring fraction, `vr` ring speed / core thermal speed, `u` electron drift /
electron thermal speed, `tau = Te/Ti`, `aniso = T⊥/T∥` applied to the *ring*
component's Maxwellian widths (used only by the moment-matched control)."""
function plasma(; delta = 1.0, vr = 0.0, u = 0.0, tau = 1.0, aniso = 1.0)
    s = NormalizedSpecies[]
    delta < 1 && push!(s, NormalizedSpecies(1.0, (1 - delta) * PI2, Maxwellian(vth = EPS)))
    if delta > 0
        wperp = EPS * sqrt(aniso)
        vdf = iszero(vr) ? Maxwellian(vth = (wperp, EPS)) :
              Maxwellian(vth = (wperp, EPS), vr = vr * EPS)
        push!(s, NormalizedSpecies(1.0, delta * PI2, vdf))
    end
    we = EPS * sqrt(tau * MI_ME)
    push!(s, NormalizedSpecies(-MI_ME, PI2 * MI_ME, Maxwellian(vth = we, vd = u * we)))
    return Tuple(s)
end

"""Bi-Maxwellian with the same n, T⊥, T∥ as the ring plasma at `(delta, vr)`.

The gyro-ring has ⟨v⊥²⟩ = w² + vr², so the moment-matched anisotropy of the heated
component is A = 1 + (vr/w)². Same first two velocity moments, different shape."""
moment_matched(; delta = 1.0, vr = 0.0, kw...) =
    plasma(; delta, vr = 0.0, aniso = 1 + vr^2, kw...)

# ---------------------------------------------------------------------------
# Growth-rate driver
# ---------------------------------------------------------------------------

const KRHO = 0.5:0.25:3.5
const RATIO = (1.0e-4, 0.005, 0.01, 0.02, 0.035, 0.05, 0.08, 0.12)   # k∥/k⊥
const BOX = (0.1 - 0.02im, 5.0 + 0.5im)                              # ω survey box, Ω_O⁺

const KGRID = [(kr, r) for kr in KRHO for r in RATIO]

wavenumber(kr, r) = Wavenumber(kr / EPS, r * kr / EPS)

"""Peak growth over the (k⊥ρ, k∥/k⊥) grid from seedless surveys.
Returns `(γ, k⊥ρ, k∥/k⊥, ω_r)` in units of Ω_O⁺."""
function gamma_max(pl; kgrid = KGRID)
    out = fill((-Inf, NaN, NaN, NaN), length(kgrid))
    Threads.@threads for i in eachindex(kgrid)
        kr, r = kgrid[i]
        sol = solve(DispersionProblem(pl, BOX, wavenumber(kr, r); mode = :P))
        best = (-Inf, NaN, NaN, NaN)
        for b in sol.roots
            w = b.omega
            isfinite(w) && real(w) > 0.15 && imag(w) > best[1] &&
                (best = (imag(w), kr, r, real(w)))
        end
        out[i] = best
    end
    return maximum(out)
end

"""Re-solve a survey root with a seeded Muller solve — a growing root that the
polisher does not return to the same place is not trusted."""
function confirm(pl, gamma, kr, r, wr; tol = 1.0e-3)
    isfinite(gamma) || return false
    w = solve(DispersionProblem(pl, complex(wr, gamma), wavenumber(kr, r); mode = :P)).omega
    return isfinite(w) && abs(w - complex(wr, gamma)) < tol * max(1.0, abs(wr))
end

"""Bisect the scalar parameter of `build` whose peak growth equals `target`.
Returns `(x_c, γ, k⊥ρ, k∥/k⊥, ω_r)`; `x_c = lo` when the plasma is already unstable there,
`NaN` when `hi` does not destabilise it."""
function threshold(build; target = 1.0e-3, lo = 0.0, hi = 0.5, iters = 9)
    g0 = gamma_max(build(lo))
    g0[1] > target && return (0.0, g0...)
    gamma_max(build(hi))[1] > target || return (NaN, NaN, NaN, NaN, NaN)
    res = (NaN, NaN, NaN, NaN)
    for _ in 1:iters
        mid = (lo + hi) / 2
        g = gamma_max(build(mid))
        g[1] > target ? (hi = mid; res = g) : (lo = mid)
    end
    return (hi, res...)
end

writetsv(path, header, rows) = open(path, "w") do io
    println(io, header)
    for r in rows
        println(io, join(r, "\t"))
    end
end
