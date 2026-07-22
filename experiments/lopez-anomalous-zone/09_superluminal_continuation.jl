# The A/IC branch beyond the light line (μ = 2), via the corrected López
# continuation of 06 followed through the gap 1 < Re z < √(1+t²).
#
# Why this is the only evaluator that can go there: for damped superluminal ω
# every real-sliced momentum integral misses the resonance-ellipse apex branch
# cut (docs/src/relativistic.typ) — VMD's closed form and CoupledVDF both warn and
# return the straight integral. AAA extrapolation from UHP samples also fails
# here (measured below): the root sits ≈3 convergence radii past the ω = k∥c
# branch point, and the fit spends its degree tracing the cut strings at
# ω = k∥c and ω = √(k∥²c²+Ω²). The corrected closed form needs no
# extrapolation: its formula is analytic along the whole root path.
include("06_corrected_continuation.jl")
using Printf

# --- holomorphy of the extension where no other evaluator is valid ---
println("\n== 09: superluminal extension ==")
println("holomorphy defect |∂z̄|/|∂z| of ΛL_corr at superluminal damped ω (k=2.5):")
for ω in (2.53 - 0.2im, 2.577 - 0.464im, 2.6 - 0.7im)
    zb, zz = dzbar(ω_ -> ΛL_corr(ω_, 2.5, 2.0), ω)
    @printf("  ω=%+.3f%+.3fim  defect=%.1e\n", real(ω), imag(ω), abs(zb) / abs(zz))
end

# --- overlap validation: subluminal k where VMD's determinant is exact ---
vmd_ref = ((1.6, 1.5482 - 0.48397im), (1.8, 1.7831 - 0.48021im))
println("\nsubluminal overlap vs VMD determinant roots:")
for (k, ωv) in vmd_ref
    r = muller(ω -> ΛL_corr(ω, k, 2.0), ωv - 1.0e-3, ωv, ωv + 1.0e-3im)
    @printf("  k=%.1f  ΛL_corr root %+.5f%+.5fim  |Δ| to VMD = %.1e\n", k, real(r), imag(r), abs(r - ωv))
end

# --- trace through and beyond the light line ---
ks = collect(1.6:0.05:3.0)
ωs = ComplexF64[]
w = 1.5482 - 0.48397im
for k in ks
    f = ω -> ΛL_corr(ω, k, 2.0)
    global w = muller(f, w - 1.0e-3, w, w + 1.0e-3im)
    push!(ωs, w)
end
println("\nA/IC continued (μ=2): crossing ωr = k∥c near k∥ ≈ 1.9, then slightly superluminal:")
for (k, ω) in zip(ks, ωs)
    k in (1.6, 1.9, 2.0, 2.2, 2.5, 3.0) &&
        @printf("  k=%.2f  ω=%+.5f%+.5fim  ωr/k=%.4f\n", k, real(ω), imag(ω), real(ω) / k)
end

open(joinpath(@__DIR__, "figdata_superluminal.csv"), "w") do io
    println(io, "k,wr,gm")
    for (k, ω) in zip(ks, ωs)
        println(io, "$k,$(real(ω)),$(imag(ω))")
    end
end

# --- far tail: damping decays slowly; no sharp cutoff at the naive resonance
# band edge ωr²−k² = Ω² (a real-ω criterion; |γ| ~ 0.4 broadens it away) ---
println("\nfar tail:")
for k in 3.2:0.2:6.0
    f = ω -> ΛL_corr(ω, k, 2.0)
    global w = muller(f, w - 1.0e-3, w, w + 1.0e-3im)
    k in (4.0, 5.0, 6.0) && @printf("  k=%.1f  ω=%+.4f%+.4fim  ωr²−k²=%.3f\n", k, real(w), imag(w), real(w)^2 - k^2)
end
println("csv written")

# --- the second damped family across the resonance band ---
# Muller needs Δk ≤ 0.02 near its light-line crossing (k ≈ 1.73) or it jumps
# onto the near-real EM branch. It traverses the whole band
# 0 < ωr²−k² < Ω² with γ ≈ −0.19 nearly flat; at the band edge (k ≈ 6.1,
# ωr²−k² = Ω²) even the optimal relativistic Doppler shift can no longer
# satisfy γ_L(ω − k∥v∥) = Ω — min over |v∥|<c is √(ω²−k∥²c²) — so the damping
# shuts off and the root lands on the real superluminal EM branch,
# ω² − k² → 2Π²K₁(μ)/K₂(μ) = 1.102 at μ = 2.
w2 = 0.92445 - 0.19622im   # VMD root at k = 0.95
println("\nsecond family across the band:")
for k in 0.96:0.02:6.0
    f = ω -> ΛL_corr(ω, k, 2.0)
    global w2 = muller(f, w2 - 1.0e-3, w2, w2 + 1.0e-3im)
    k in (1.0, 1.73, 2.0, 3.0, 4.0, 5.0, 6.0) &&
        @printf("  k=%.2f  ω=%+.5f%+.5fim  ωr²−k²=%.4f\n", k, real(w2), imag(w2), real(w2)^2 - k^2)
end
for k in (6.4, 7.0)
    f = ω -> ΛL_corr(ω, k, 2.0)
    global w2 = muller(f, w2 - 1.0e-3, w2, w2 + 1.0e-3im)
    @printf("  k=%.2f  ω=%+.5f%+.5fim  ωr²−k²=%.4f  (post-exit: on the EM branch)\n", k, real(w2), imag(w2), real(w2)^2 - k^2)
end

# --- small-k end: the aperiodic ladder ---
# ΛL is a single branch (no R/L squaring), real on the imaginary axis, so its
# axis roots are plain sign changes — no double-zero obstruction. The "second
# family" is the shallow end of a ladder that densifies as k → 0 (the k = 0
# relativistic cyclotron continuum ω = Ω/γ_L approached by discrete damped
# roots); the γ = −1.27 member is the deep family of the Fig. 5 replica.
function ladder(k; smax=1.5, n=3000)
    f = s -> real(ΛL_corr(complex(0.0, -s), k, 2.0))
    ss = range(0.002, smax, length=n)
    vals = map(f, ss)
    out = Float64[]
    for i in 1:(n-1)
        sign(vals[i]) == sign(vals[i+1]) && continue
        lo, hi = ss[i], ss[i+1]
        for _ in 1:50
            m = (lo + hi) / 2
            f(m) * vals[i] > 0 ? (lo = m) : (hi = m)
        end
        push!(out, -(lo + hi) / 2)
    end
    return out
end
println("\naperiodic ladder (axis roots of ΛL_corr):")
for k in (0.02, 0.05, 0.1, 0.2, 0.3)
    println("  k=", k, "  γ: ", round.(ladder(k), digits=4))
end

# --- the k → 0 ladder is exact: γ_n = −2μΩ/(π(2n−1)) ---
# (numerically to ~1e-7 at k = 0.001, for μ = 2…25 and n = 1…3; the deep
# "aperiodic family" of the Fig. 5 replica is the n = 1 member.)
function axis_root(μ, guess; k=0.001)
    f = s -> real(ΛL_corr(complex(0.0, -s), k, μ))
    lo, hi = 0.9guess, 1.1guess
    for _ in 1:60
        m = (lo + hi) / 2
        f(m) * f(lo) > 0 ? (lo = m) : (hi = m)
    end
    return -(lo + hi) / 2
end
println("\nk→0 ladder law γ_n = −2μΩ/(π(2n−1)):")
for μ in (2.0, 10.0, 25.0), n in (1, 2)
    pred = 2μ / (π * (2n - 1))
    got = axis_root(μ, pred)
    @printf("  μ=%-4.0f n=%d  measured %+.5f  law %+.5f  Δ=%.1e\n", μ, n, got, -pred, abs(got + pred))
end

# --- derivation check (report §5.1): as y → 0 on the axis,
#   ΛL_corr(−is, y) → D(s) − (2πs/(μK₂y³))·cos(μ/s)
# The y⁻³ term is the σ=2 continuation tail with its endpoints at the
# continued resonant energies γ_res → ∓1/x = ±i/s (elementary Gaussian-free
# integrals give −(4x²/μ³)cosh(μ/x) → cos on the axis); D(s) is the O(1)
# θ-less background after the O(t⁻¹) moment cancels the −μ/y² cold term.
# Zeros are pinned to cos(μ/s) = 0 ⇒ the ladder.
function Dbg(s, μ)
    I, _ = quadgk(γ -> exp(-μ * γ) * γ * (γ^2 - 1)^1.5 / (s^2 * γ^2 + 1), 1.0, Inf; rtol=1.0e-10)
    return 1 - μ^2 / (3besselk(2, μ)) * I
end
Ccont(s, μ, y) = -(2π * s / (μ * besselk(2, μ) * y^3)) * cos(μ / s)
println("\nasymptotic form ΛL = D(s) + C(s,y):")
for (μ, s) in ((2.0, 1.0), (10.0, 2.5)), y in (0.02, 0.01, 0.005)
    lhs = real(ΛL_corr(complex(0.0, -s), y, μ))
    rhs = Dbg(s, μ) + Ccont(s, μ, y)
    @printf("  μ=%-4.0f s=%.1f y=%.3f  rel.err=%.1e (→ 0 like y)\n", μ, s, y, abs(lhs - rhs) / abs(lhs))
end
# López's frozen-support θ term instead starts at γ₁(Re z=0) ≈ 1/y: their
# continuation is O(e^{−μ/y}) — the whole aperiodic ladder is absent from the
# non-holomorphic closed form (their Λ stays smooth and positive through
# s = 4/π at k = 0.01 while the true function swings ±3e6 through zero).

# --- finite-k continuation of the ladder: the in-band quasimode stack ---
# Ladder members leave the axis pairwise and become propagating damped modes
# hugging the light line; member count at fixed k scales with μ like the
# ladder. (μ=2's "second family" is the least-damped member.)
println("\nquasimode stack at k = 1.5 (Muller seed grid, |Λ| < 1e-9):")
for μ in (2.0, 10.0)
    found = ComplexF64[]
    for wr in 0.1:0.2:1.7, wi in (-0.3, -0.7, -1.1, -1.6)
        r = muller(ω -> ΛL_corr(ω, 1.5, μ), complex(wr, wi) - 1.0e-3, complex(wr, wi), complex(wr, wi) + 1.0e-3im)
        (isfinite(r) && abs(ΛL_corr(r, 1.5, μ)) < 1.0e-9 && 0.01 < real(r) < 2.5 && -2.5 < imag(r) < -0.01) || continue
        any(z -> abs(z - r) < 1.0e-3, found) || push!(found, r)
    end
    println("  μ=", μ, "  (incl. the A/IC root): ", [round(z, sigdigits=4) for z in sort(found, by=imag, rev=true)])
end
