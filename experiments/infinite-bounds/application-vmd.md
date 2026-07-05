# Application: bounds-free CoupledVDF (PR #36, open)

How the sinc/cot method (README.md) is wired into VlasovMaxwellDispersion.

The per-harmonic parallel integrals of the susceptibility are exactly the
README's problem with g = the 5-vector of gradient slice moments
`(∂⊥f, u∂⊥f, u²∂⊥f, ∂∥f, u∂∥f)` at fixed p⊥, and ζ_n = (ω − nΩ)/k∥ a ladder.
The finite `para`/`perp` box the API used to require was an artifact of the
old log-based Plemelj subtraction (divergent on ℝ), not physics. PR #36 makes
bounds optional (default p∥ ∈ ℝ, p⊥ ∈ (0,∞)):

- `src/distributions/CoupledVDF.jl`: `_coupled_perp_sinc` = the ladder
  evaluator fused with the perp Bessel bilinears; boundedness is a TYPE
  parameter on `CoupledVDF`, so finite boxes keep the old adaptive path as a
  separate specialization (no TTFX coupling). Constructor hoists the f₀-only
  moments (⟨p⊥²⟩ → harmonic window; ⟨p∥⟩, ⟨p∥²⟩ → map center/scale u_c, S).
- Results: 3–6× faster than a hand-tuned box at ≤1e-7 vs analytic oracles
  (Maxwellian k⊥=1: 1.3 vs 4.5 ms; BiKappa κ=2: 1.7 vs 8.9 ms — where a ±30
  box silently truncates to 2.4e-5; ring vr/vth⊥=12: 47 vs 92 ms). Cold
  compile guard ~0.8× main.
- Finite bounds remain required for `regime=Relativistic()` (fixed-GL box),
  `closure=Newberger()`, and GridVDF splines (not strip-analytic).
- `GyroRing` rewritten to the scaled-Bessel form
  `exp(-((v−vr)/vth)²)·besselix(0,·)` — raw `besseli` overflows once infinite
  tails probe large v (README §2 overflow rule).

## Traps hit during integration (do not rediscover)

- **Zero-integral quadgk hang**: ⟨p∥⟩ is exactly 0 for symmetric f₀; a scalar
  adaptive quadrature can never meet a relative tolerance on a zero integral
  and refines to maxevals ("constructor stuck") — the constructor computes
  the moment family as ONE fused SVector quadrature so the vector norm lends
  the zero component a scale. Keep it fused.
- **quadgk(-Inf,-S,S,Inf) breakpoint seeding** crashed the Julia session
  (undiagnosed; guard with maxevals before trusting).
- **Mapping/scaling the outer perp (0,∞) integral** regressed (ring 243 vs
  164 ms) — QuadGK's builtin transform + adaptivity beats a global stretch.
  The perp direction stays adaptive.

## Open threads

- The same cot identity applies per-p⊥-slice to the relativistic rationalized
  resonance poles (docs/relativistic.md) — could replace fixed-GL peeling;
  the γ-artifact roots and the apex branch cut need thought.
- SeparableVDF / ReducedVDF's 1-D parallel `hilbert` (integrals.jl) is the
  same scalar problem — sinc/cot could replace adaptive+subtraction there.
- `infbounds.jl` keeps the VMD-coupled prototypes of the superseded
  approaches (README §6) plus the sinc prototype, runnable against package
  internals.
