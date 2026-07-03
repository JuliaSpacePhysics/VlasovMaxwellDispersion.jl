# VlasovMaxwellDispersion.jl — Architecture

Arbitrary gyrotropic VDFs (analytic function forms or numerical grids). Dimensionless internally.

## 1. The reduction (the spine)

For any gyrotropic `f₀` the susceptibility starts from one helical-orbit
phase-space integral. The default evaluator takes the textbook harmonic
expansion; then every term, when `f₀` is separable, factors into a
**parallel** Cauchy/Hilbert integral `H∥` (pole at `ζ_n=(ω−nΩ)/k∥`) times a
**perpendicular** Bessel moment `P⊥`. Two 1-D primitives carry the separable
core; the only variation is how their integrand is represented:

| representation | `H∥` (parallel) | `P⊥` (perp) |
|---|---|---|
| Maxwellian (Gaussian) | `√π Z(ζ)` | `Γ_n(λ)=Iₙ(λ)e^{−λ}` |
| piecewise polynomial | log-ratio formula (branch-cut invariant) | Bessel-product power series (Schläfli ₂F₃) per cell |

Maxwellian/Kappa are fast overloads.

### 1a. Grid → basis projection (the other hard half)

Turning a tabulated/noisy `f₀` into the piecewise-poly coefficients above is its
own step and must guarantee what the spine assumes: **f₀ ≥ 0**, **C¹** (so
∂f₀/∂v feeding the residue is smooth), and a **controlled approximation error**.
Default: **non-negative least-squares tensor B-spline** — adaptive knots to a
relative-error tolerance, positivity by construction, analytic derivatives. Naïve
per-cell interpolation gives negative f₀ and noisy ∂f₀ near the pole and is
rejected.

## 2. Traits (orthogonal dispatch)

| trait | values | drives |
|---|---|---|
| `Regime` | `NonRelativistic` / `Relativistic` | active coordinates, `γ`, pole map |
| `IntegralClosure` | `HarmonicSum` / `Newberger` |

Specializations are trait combinations: 
- Maxwellian/Cold = `Analytic+Separable+NonRel+HarmonicSum`
- Maxwell–Jüttner = `Analytic+Relativistic`
- gridded relativistic `f₀(p̂⊥,p̂∥)` = `PiecewisePoly+Coupled+Relativistic`.

## 3. Cross-validation

Three reference solvers in `external/` serve as ground truth for tests:

| reference | validates |
|---|---|
| `LinearMaxwellVlasov.jl` | bi-Maxwellian χ numbers, Newberger coupled path, complex-k |
| `ALPS` | arbitrary gyrotropic + **relativistic** test inputs |
| `MPDES` | piecewise-poly `H∥`/`P⊥`, NNLS spline, GES global finder, paper figures |

Plus analytic anchors with no external dep: Stix cold R/L/O/X, Maxwellian→cold limit, Langmuir+Landau vs the `Z`-function dispersion, electrostatic limit.

accelerations, both learned from the references:

1. *Velocity integral* — replace nested adaptive QuadGK with precompute-once:
   either MPDES-style **project `f₀`→2-D piecewise-poly then closed-form per cell**
   (`projection.jl`+`hilbert_pwpoly.jl` all exist; analytic in ω,
   AD-clean) or ALPS-style **fixed-grid Simpson + precomputed Bessel weights**.
   This fixes the cost that dominates even at small `k⊥ρ`.

## 4. Implementation status

**Built + validated**:

All kinetic paths specialize one nested 2-D integral (`derivation.md`):
`CoupledVDF` (general) ⊃ `SeparableVDF` (factors) ⊃ Maxwellian/MJ (`Z`/`Γ_n`
closed form).

| capability | file | validation |
|---|---|---|
| `H∥` piecewise-poly + branch cut | `hilbert_pwpoly.jl` | exactness + Plemelj (17) |
| `P⊥` Bessel×poly | `perp_analytic.jl` | vs quadrature (11) |
| NNLS B-spline projection | `projection.jl` | positivity/C¹/error-bound |
| relativistic Maxwell–Jüttner | `MaxwellJuttner.jl` | ALPS `relativistic` roots 1&2; μ→∞ limit |
| GRPF `GlobalDispersionProblem` + `BranchProblem` | `solve.jl`/`track.jl` | ALPS `kpar_fast` scan (rtol 1e-2); GRPF unit |