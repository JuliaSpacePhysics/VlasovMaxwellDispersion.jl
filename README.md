# VlasovMaxwellDispersion.jl

Linear hot-magnetized Vlasov–Maxwell dispersion solver. 

Full derivation: [`derivation.md`](docs/derivation.md).
Method and internals: [`architecture.md`](docs/architecture.md).

## Quick start

```julia
using Pkg; Pkg.develop(path="."); Pkg.instantiate()
using VlasovMaxwellDispersion

pl = Plasma(Species(-1.0, 1.0, Maxwellian(1.0)))   # Ω̃, Π̃², VDF
k  = Wavenumber(0.0, 0.7)                           # k̃ = (k⊥, k∥)·c/Ω_ref

alg = Muller()
sol = solve(LocalDispersionProblem(pl, k, 1.2 - 0.1im), alg)        # local root (Langmuir+Landau)
ω   = sol.omega

gsol = solve(GlobalDispersionProblem(pl, k, (0.5 - 0.6im, 2.5 + 0.1im)))  # all roots+poles in a box
roots, poles = gsol.omega, gsol.poles

ks   = [Wavenumber(0.0, kz) for kz in 0.3:0.05:1.0]
ωs   = solve(BranchProblem(pl, ks, 1.2 - 0.1im)).omega             # track a branch in k
```

Arbitrary distribution — e.g. bump-on-tail instability:

```julia
f(u) = 0.94exp(-u^2)/√π + 0.06exp(-((u-4)/0.5)^2)/(0.5√π)
pl = Plasma(Species(-1.0, 1.0, SeparableVDF(f; lower=-12.0, upper=14.0)))
ω  = solve(LocalDispersionProblem(pl, Wavenumber(0.0, 0.25), 1.0 + 0.05im)).omega   # Im ω > 0 ⇒ growth
```

## Capabilities

| VDF / mode | path |
|---|---| 
| Arbitrary Analytic `f₀(p∥,[⊥)` |
| Grid / numerical VDF | NNLS B-spline → piecewise-poly `H∥`,`P⊥` |
| Maxwell–Jüttner (relativistic) | Trubnikov/Swanson integral |
| (bi-)Maxwellian, drifting | `Z`-function harmonic sum |
| Cold fluid | Stix S,D,P closed form |

Solvers follow the `CommonSolve.solve(problem, algorithm)` interface: a
`Local`/`Global`/`BranchProblem` (seed / search box / k-sequence) solved by
`Muller` (default; `Secant` alternative) / `GRPF` (`RootsAndPoles.jl`) /
`ArcLength`. `solve` returns a `DispersionSolution` (`.omega`, `.poles`,
`.retcode`).

Two closure for orbit integrals are available (`derivation.md` §3):

- `HarmonicSum()` (default): truncate harmonic sum with `nmax ≈ k⊥ρ`, one Landau `hilbert` per harmonic. Handles damping via the Landau contour.
- `Newberger()` — use Qin closed-orbit `T(a,z)` in complex-order Bessel `J_{±a}`, with no truncation (cost flat in `k⊥ρ`), Handles damping with residue extraction.