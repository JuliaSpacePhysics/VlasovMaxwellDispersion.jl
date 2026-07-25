# VlasovMaxwellDispersion.jl — Architecture

Arbitrary gyrotropic VDFs (analytic function forms or numerical grids). Dimensionless internally.

Susceptibility starts from one helical-orbit phase-space integral.
The default evaluator takes the textbook harmonic expansion.

## Traits (orthogonal dispatch)

All kinetic paths specialize one nested 2-D integral (`derivation.md`):
`CoupledVDF` (general) ⊃ `SeparableVDF` (factors) ⊃ Maxwellian/MJ (`Z`/`Γ_n`
closed form).

| trait | values | drives |
|---|---|---|
| `Regime` | `NonRelativistic` / `Relativistic` | active coordinates, `γ`, pole map |
| `IntegralClosure` | `HarmonicSum` / `Newberger` | harmonic truncation vs. closed-orbit `T(a,z)`; damping via Landau contour vs. residue extraction |
| `ChiBackend` | `FixedNodeEval` / `AdaptiveEval` | whether a per-`k` plan hoists the ω-independent work onto fixed nodes, or the adaptive quadrature re-runs at every ω |

Specializations are trait combinations: 
- Maxwellian/Cold = `Analytic+Separable+NonRel+HarmonicSum`
- Maxwell–Jüttner = `Analytic+Relativistic`
- gridded relativistic `f₀(p̂⊥,p̂∥)` = `PiecewisePoly+Coupled+Relativistic`.

## Cross-validation

Three reference solvers in `external/` serve as ground truth for tests:

| reference | validates |
|---|---|
| `LinearMaxwellVlasov.jl` | bi-Maxwellian χ numbers, Newberger coupled path, complex-k |
| `ALPS` | arbitrary gyrotropic + **relativistic** test inputs |
| `MPDES` | piecewise-poly `H∥`/`P⊥`, NNLS spline, GES global finder, paper figures |

Plus analytic anchors with no external dep: Stix cold R/L/O/X, Maxwellian→cold limit, Langmuir+Landau vs the `Z`-function dispersion, electrostatic limit.

## The per-`k` plan layer

The nested adaptive QuadGK over `(p⊥, p∥)` dominates even at small `k⊥ρ`, and at fixed
`k` the ONLY ω-dependence is the Landau pole `ζₙ=(ω−nΩ)/k∥`. So every kinetic path splits
into a `plan_contribution(species, k)` that hoists the ω-independent work, and a
`plan(ω)` that is cheap arithmetic. Two accelerations from the references, both wired: