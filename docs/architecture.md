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

Acceleration learned from the references:

1. *Velocity integral* — replace nested adaptive QuadGK with precompute-once:
   either MPDES-style **project `f₀`→2-D piecewise-poly then closed-form per cell**
   (`projection.jl`+`hilbert_pwpoly.jl` all exist; analytic in ω,
   AD-clean) or ALPS-style **fixed-grid Simpson + precomputed Bessel weights**.
   This fixes the cost that dominates even at small `k⊥ρ`.
