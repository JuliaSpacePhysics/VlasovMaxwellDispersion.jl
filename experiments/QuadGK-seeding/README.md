# CoupledVDF QuadGK seeding: compile-latency regression

Investigation of a seeding experiment which seeds the nonrelativistic harmonic quadratures 
with initial breakpoints to cut adaptive over-refinement of f₀'s near-zero tails.

Outcome: its cold-compile cost (~4×) outweighed the warm win for the suite/CI. 

The change itself is preserved as `experiments/QuadGK-seeding/coupled_seeding.patch`.

### Cold latency (`@timed(...).compile_time`, fresh process, s) — regressed

Single coupled `contribution` call: **OLD 3.4s → NEW 14.6s**.
Batch over all nonrel coupled cases (TTFX, fresh process): **11.9s → 38.3s (+26s)**.

## Root cause

QuadGK's initial-segment handling is unrolled over the breakpoint tuple, so
**compile time grows with the number of initial breakpoints.** Two seeds were
added; ablation by swapping the source and remeasuring batch TTFX:

| variant | TTFX | Δ vs OLD |
|---|---|---|
| OLD — no seed | 11.9s | — |
| inner `_segs` Val(6) only (plain outer) | 18.4s | +6.5s |
| + outer `_osegs` variable = **NEW (shipped)** | 38.3s | +26s |
| inner Val(6) + outer **static** Val(8) | 44.0s | +32s |
| inner Val(3) + outer static Val(4) | 31.9s | +20s |

Takeaways:
- The **outer perp seed dominates** (~+20s). The heavy `_coupled_perp` closure is
  recompiled for the seeded outer `quadgk`, and it is costly per breakpoint.
- It is *not* mainly the variable arity: a single static Val(8) outer seed is
  *worse* than the variable one. **Breakpoint count is the lever**, not
  static-vs-variable. (Variable arity does add its own per-distinct-count
  specializations on top.)
- The inner static seed is comparatively cheap (+6.5s) and carries most of the
  warm-runtime win.

## Mitigation options (decision deferred)

| option | cold Δ | warm gain kept | cost |
|---|---|---|---|
| A. Drop outer `_osegs`, keep inner `_segs` | +26s → +6.5s | ~70% | trivial, no deps |
| B. Drop both seeds (revert) | → 0 | 0 | loses the commit |
| C. Keep both, make outer static + PrecompileTools workload | → ~0 at test time | 100% | +dep, +workload; compile moves into cached `buildpkg` |

A is the pragmatic default. C keeps full warm speed but needs a static outer
arity first (PrecompileTools cannot bake an unbounded set of `_osegs` arities)
and a precompile workload that exercises the coupled path.

## Reproduce

`compile_time` of a first coupled call in a fresh process is the metric. The probe
script measures whatever is checked out:

```
julia --project=. experiments/QuadGK-seeding/coupled_compile_latency.jl     # main: ~3–4 s/case
git apply experiments/QuadGK-seeding/coupled_seeding.patch                  # re-add the seeding
julia --project=. experiments/QuadGK-seeding/coupled_compile_latency.jl     # seeded: ~15–17 s/case
git checkout -- src/distributions/CoupledVDF.jl                  # undo
```

Related test:

```julia
using VlasovMaxwellDispersion, TestItemRunner
@run_package_tests filter = ti -> (:latency in ti.tags)
```
