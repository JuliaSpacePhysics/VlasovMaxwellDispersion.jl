# Sheared E×B flow: what reduces to an eigenvalue problem, and when it is even kinetic

Single source of truth for this experiment. Code files carry no derivations, only pointers
to the equation numbers here.

Setup throughout: `B₀ = B ẑ`, static transverse field `E₀(x) x̂`, and the E×B drift it
implies, `u(x) = (E₀×B)/B²|_y`. Species have signed gyrofrequency `Ω = qB/m`. The one
control parameter is

  **σ ≡ S/Ω,  S ≡ du/dx**  (species-signed: σₑ ≈ (mₑ/mᵢ)σᵢ ≈ 0, so shear is an ion effect)

Electrostatic throughout — EIC, EIH and IEDDI all live there — but §1–§4 are statements
about orbits, so they carry to the full tensor unchanged.

---

## §0. Is this a kinetic problem at all?

Ask before computing. Rewrite σ using `S ~ Δu/L`, `ρ = w/Ω`, `w = √(2T/m)`:

  **σ = (Δu/w)·(ρ/L)**   (1)

That identity governs everything. It couples "how strong is the shear" to "how narrow is
the layer", and the two cannot be separated.

| regime | ordering | fastest mode | right tool |
|---|---|---|---|
| broad layer, subsonic | `L ≫ ρ`, `Δu ≲ w` ⇒ `σ ≪ 1` | Kelvin–Helmholtz, `γ ≈ 0.2 S`, `ω ≈ kᵧΔu ≪ Ω`, `kᵧρ ≪ 1` | **fluid.** MHD/two-fluid is correct and cheaper |
| gyroharmonic band | any σ, `ω ≈ nΩᵢ` | EIC / IEDDI, ion-cyclotron harmonics | **kinetic, mandatory.** No fluid model has this branch |
| narrow layer | `L ≲ few ρᵢ` ⇒ `σ = O(1)` | EIH, `γ ~ Ω`, `kᵧρ ~ 1` | **kinetic and nonlocal** (§5, §6) |

Two consequences worth stating plainly.

**(i) For the ionospheric target the fluid question does not arise.** The observed auroral
waves are O⁺ cyclotron harmonics at tens of Hz, `ω ≈ nΩ_O⁺`. A fluid model has no
gyroharmonic branch to destabilise — not a quantitative error, an absent mode. Whatever
the shear does there, only a kinetic calculation can see it.

**(ii) Strong shear and locality are incompatible.** By (1), `σ = O(1)` with `Δu` limited to
a few `w` — which the ionosphere is, `Δu ~ 1–5 km/s` against `w(O⁺, 0.3 eV) = 1.9 km/s` —
forces `L ~ ρᵢ`. So there is no regime that is simultaneously strongly sheared and locally
describable. §4 measures exactly this, and the measured error is what (1) predicts.

Auroral numbers, `B = 2.8×10⁻⁵ T`, `T(O⁺) = 0.3 eV`: `ρ_O⁺ = 11 m`, `Ω_O⁺/2π = 27 Hz`.
Layer widths `L ~ 0.1–1 km` and `Δu ~ 1–5 km/s` give **σ ≈ 0.01–0.3** — straddling the
boundary. Fluid KH (`γ ≈ 0.2S`) overtakes a typical EIC growth (`γ ~ 0.05 Ω`) only above
`σ ≈ 0.25`, i.e. only for the narrowest layers, which are also the ones where a local
kinetic treatment has already failed.

---

## §1. Orbits, invariants, circularizing coordinates

Take **linear** shear, `u(x) = u₀ + Sx`, `S` constant. Equations of motion:

  `v̇ₓ = Ω(v_y − u(x)),  v̇_y = −Ω vₓ,  v̇_z = 0`

`v̇_y = −Ω vₓ` integrates to a conserved `P = v_y + Ωx`; eliminating `v_y` gives a harmonic
oscillator,

  `ẍ = −Ω²(1+σ)(x − X),   X ≡ (v_y + Ωx − u₀)/(Ω(1+σ))`   (2)

so orbits stay closed and harmonic at the **effective gyrofrequency**

  `Ω_eff ≡ Ω√(1+σ)`   (3)

with `X` the conserved guiding centre. Eliminating back,

  `v_y − u(X) = −Ω(x − X)`   (4)

— guiding centres drift at exactly the local E×B speed, and the gyration is an **ellipse**:
`vₓ` and `v_y − u(X)` oscillate in quadrature at `Ω_eff` with amplitude ratio `(1+σ)^(−1/2)`.
The second perpendicular invariant is

  `w² ≡ vₓ² + (1+σ)(v_y − u(X))²`   (5)

Invariants `(w², X, v_z)` ⇒ any `f₀ = F(w², X, v_z)` is an exact Vlasov equilibrium.
`∂F/∂X = 0` is pure velocity shear (uniform n, T); `∂F/∂X ≠ 0` is how density and
temperature gradients would enter. Everything below takes `∂F/∂X = 0`.

Put

  `ξ ≡ vₓ,  η ≡ √(1+σ)(v_y − u(X)),  w² = ξ² + η²`   (6)

By (4) the orbit in `(ξ,η)` is a **circle** of radius `w` traversed at `Ω_eff`, same sense as
`Ω`. The shear is a linear map on velocity space.

**Corollary (forced anisotropy).** At fixed `x`, (5) reads `w² = vₓ² + (v_y−u₀)²/(1+σ)`: the
lab-frame variance along the flow exceeds that across it by `1+σ`. No isotropic sheared
Vlasov equilibrium exists; this is forced, not assumed.

## §2. The step that is easy to get wrong

`X` depends on `v_y`, so `∂F/∂v_y` is not the naive chain rule on (5). With
`∂X/∂v_y = 1/(Ω(1+σ))`:

  `∂w²/∂v_y = 2(1+σ)(v_y−u)·(1 − S ∂X/∂v_y) = 2(1+σ)(v_y−u)·(1 − σ/(1+σ)) = 2(v_y−u)`

The `(1+σ)` cancels, leaving

  `∂F/∂vₓ = 2F′ξ,   ∂F/∂v_y = 2F′η/√(1+σ)`   (7)

Miss the guiding-centre dependence and the numerator picks up `√(1+σ)` instead of its
inverse; numerator and orbit phase then carry different effective wavevectors and the
problem looks irreducible. It is not — they agree, and that is what makes §3 work.

## §3. Exact reduction of the orbit part

Perpendicular displacement along the orbit (`r(0)=0`):

  `Δx(τ) = X(1−cos Ω_eff τ) + (vₓ/Ω_eff) sin Ω_eff τ`
  `Δy(τ) = u(X)τ + (ΩX/Ω_eff) sin Ω_eff τ − (Ωvₓ/Ω_eff²)(1−cos Ω_eff τ)`   (8)

The oscillatory part of `k·Δr` has amplitude `(w/Ω_eff)·k_eff`, and by (7) the numerator
`k·∂F/∂v|_⊥ = 2F′(κ·(ξ,η))` carries the **same** vector:

  `κ = (kₓ, k_y/√(1+σ)),   k_eff ≡ |κ| = √(kₓ² + k_y²/(1+σ))`   (9)

The velocity-space Jacobian `d³v = dξ dη dv_z/√(1+σ)` cancels against the renormalisation
of `F`, so `Π` is untouched. Dropping for a moment the secular `u(X)τ` term in (8):

  **χ_s^shear(ω,k;S) = [(k_eff²+k_z²)/k²] · χ_s^std(ω, (k_eff, k_z)) with Ω_s → Ω_s√(1+σ_s)**  (10)

`χ_s^std` is the ordinary gyrotropic susceptibility this package already computes, reading
the same VDF shape in `w`. The prefactor is the only bookkeeping: the response bracket
lives at `k_eff`, Poisson divides by the true `k`.

**Verified numerically.** `00_validate_transform.jl` integrates the response directly along
exact sheared orbits in lab coordinates — no invariants, no map, no Bessel expansion — and
compares. With the `u(X)τ` term switched off, (10) is exact to **≤1.6×10⁻⁸** (typically
10⁻¹⁰) over `σ ∈ [−0.4, 0.8]`. The reduction is a theorem, not an approximation.

So (10) is free: `Ω → Ω√(1+σ)`, `k⊥ → k_eff`, one scalar prefactor. It says the **whole
gyroharmonic ladder slides to `nΩ√(1+σ)`**, dragging the EIC resonance with it. That is a
real, purely kinetic effect with no fluid counterpart.

**But it carries no free energy.** After (6) the reduced problem is an ordinary gyrotropic
one whose equilibrium is `F(w²)` with `F′ < 0` — monotone decreasing in the invariant. No
positive slope, no inversion, nothing to feed a wave. Shear cannot destabilise anything
through (10) alone; it can only rescale what the unsheared plasma already does.

`01_eic_shear_shift.jl` confirms this on the topside auroral EIC branch (O⁺/e⁻,
`k⊥ρ = 1.5`, `k∥/k⊥ = 0.05`): the ladder tracks `√(1+σ)` to 0.2%, while the current
threshold barely moves and moves the *wrong* way for both signs of σ.

| σ | √(1+σ) | ω_r/Ω | γ/Ω at u=0.2 | u_c/v_e | u_c/u_c(0) |
|---|---|---|---|---|---|
| −0.10 | 0.9487 | 1.1237 | +0.0169 | 0.1001 | 1.009 |
| −0.02 | 0.9899 | 1.1731 | +0.0191 | 0.0991 | 0.999 |
| 0 | 1.0000 | 1.1849 | +0.0195 | 0.0992 | 1.000 |
| +0.02 | 1.0100 | 1.1963 | +0.0198 | 0.0995 | 1.003 |
| +0.10 | 1.0488 | 1.2397 | +0.0209 | 0.1012 | 1.020 |

A 5% slide of the ladder buys a 2% *rise* in threshold. **Every shear-driven instability —
KH, EIH, IEDDI, shear-modified EIC — therefore lives entirely in the term of §4.** This is
the single most consequential result here, and it inverts the implementation priority: the
cheap wrapper is a rescaling tool, not an answer to the science question.

## §4. The term that does not reduce

Restore `u(X)τ`. Since `X = x + η/(Ω√(1+σ))` is an invariant, `u(X)` is constant along an
orbit but varies across the distribution: at the analysis point `x = 0`,

  `k_y u(X) = k_y u₀ + k_y σ η₀/√(1+σ)`   (11)

so the resonance condition becomes

  `ω − k_y u₀ − (k_yσ/√(1+σ)) η₀ − k_z v_z − n Ω_eff = 0`   (12)

a **Landau-type resonance in the perpendicular invariant η₀** — guiding centres free-stream
along `y` at a speed set by where they sit. This cannot be folded into (10): `η₀ = w cos φ`
depends on the gyrophase at the observation time, so the `φ` integral no longer forces
`n = m` and the double harmonic sum does not collapse. That collapse is the entire reason
the standard kernel is cheap.

Measured cost of dropping it (`|ε_exact − ε_eq10|/|ε_eq10|`, Maxwellian ion, `ω = 1.2+0.3i`):

| σ \ kᵧρ | 0.4 | 0.8 | 1.2 | 1.6 |
|---|---|---|---|---|
| 0.05 | 0.062 | 0.073 | 0.060 | 0.040 |
| 0.10 | 0.127 | 0.155 | 0.130 | 0.090 |
| 0.20 | 0.257 | 0.309 | 0.256 | 0.179 |
| 0.40 | 0.461 | 0.483 | 0.365 | 0.253 |

Empirically **err ≈ 1.2σ**, nearly flat in `kᵧρ` across the band that matters. This is §0(ii)
in numbers: the local approximation degrades linearly in the same parameter that measures
how narrow the layer is. Caveat: this is the error on `ε_L`, not on a growth rate; a root
near marginal stability can amplify it.

Physically, (12) *is* the fluid advection term — the `n = 0` branch of (12) is the KH
critical layer. So the split is clean: **(10) is the kinetic content of shear (orbit
deformation, gyroharmonic ladder), (12) is the fluid content (flow advection), and their
interference at `n ≥ 1` is the shear-driven cyclotron-harmonic drive.** Any treatment that
keeps only one of them is describing only one channel.

## §5. Is a plane wave even a mode?

For unbounded linear shear, no. The operator contains `Sx ∂_y`, which in Fourier space is
`−iSk_y ∂_{k_x}` — first order in `k_x`. The exact solutions are Kelvin shearing modes with

  `k_x(t) = k_x(0) − S k_y t`

so wavevectors drift and the spectrum is continuous: uniform shear supports transient
amplification, not eigenvalues. A plane-wave dispersion relation is meaningful only while
`S k_y/k_x ≪ γ`, i.e. for modes fast compared with the shear itself. Combined with §4, the
honest scope of (10) is:

  **σ ≪ 1, and `γ ≫ S k_y/k_x`** — weak shear as a perturbation of a kinetic mode.

That is exactly the useful regime for the ionospheric question: `σ ≈ 0.01–0.1` (layers
`L ≳ 300 m`), where shear shifts the gyroharmonic ladder by `O(σ/2)` and moves the EIC
current threshold, with (10) accurate to `≈1.2σ`.

## §6. Narrow layers: still an eigenvalue problem, now a matrix one

For `u(x)` structured on `ρᵢ`, keep `x` in real space and Fourier only `(y,z,t)`:

  `−φ″(x) + (k_y²+k_z²)φ(x) = Σ_s (1/ε₀) ∫dx′ K_s(x,x′;ω,k_y,k_z) φ(x′)`   (13)

`K_s` is the guiding-centre response: local in `X`, smeared over the orbit, so its range is
`~2ρ_s` — banded in units of `ρᵢ`, not dense. Discretising `x` on `N` points,

  `M(ω)φ = 0,  M(ω) = D₂ + (k_y²+k_z²)I − Σ_s K_s(ω)`   (14)

`det M(ω)` is a scalar analytic function of `ω` — precisely what this package's root layer
already consumes (`Muller`, `GRPF`, `AAA` box surveys). Nothing in those solvers knows the
scalar came from a 3×3 `det 𝒟`. Building `K_s` is the work: same orbit integral, with the
`x`, `x′` endpoints no longer collapsed by a plane-wave phase. Its WKB limit must reproduce
(10) at the local `S(x)` — the cheap correctness check for any implementation.

**PIC is never needed for linear stability.** Uniform shear is algebraic (§3) within its
scope (§5); narrow layers are a banded matrix whose determinant feeds existing root
finders (§6). PIC earns its cost only at saturation.

## §7. What this implies for implementation

1. **(10) as a thin wrapper** — done, `shear.jl`, ~25 lines, no new kernels, exact to
   10⁻¹⁰ against direct orbit integration. Use it to rescale the gyroharmonic ladder and to
   check any of the below in the `σ → 0` limit. It answers no stability question by itself.
2. **Exact local `ε_L` keeping (12)** — the mandatory next step, since §3 shows the drive is
   entirely there. The direct orbit quadrature in `00_validate_transform.jl` is already the
   correct reference at ~2 s per `(ω,k)`; the fast form keeps the standard harmonic sum but
   wraps it in a gyrophase quadrature with `ω → ω − k_yu(X)` per node — one `Z` evaluation
   per `(w, φ, n)`, milliseconds, and the `n≠m` coupling that (10) has to drop is retained.
   This is what can actually answer "does shear lower the EIC current threshold".
3. **(14)** — the real generalisation, and the only route to `L ~ ρᵢ`, where §0 says the
   ionospheric shear layers actually sit. Banded `K_s`, then reuse the existing root finders
   unchanged. Its WKB limit must reproduce 2, which must reproduce 1 as `σ → 0`.

The chain 1 → 2 → 3 is also the test chain: each stage is the previous one's arbiter.
