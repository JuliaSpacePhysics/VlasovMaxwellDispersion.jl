# Development notes

[![Coverage](https://codecov.io/gh/JuliaSpacePhysics/VlasovMaxwellDispersion.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSpacePhysics/VlasovMaxwellDispersion.jl)

## Roadmap

- Optional **AAA-rational** backend for smooth analytic input (`BaryRational.jl`) with Landau-causal pole filtering.
- Check and potentially reuse `BSplineKit.jl`/`Dierckx.jl`

 1. Expose Tuned GRPF Params
  Current wrapper hardcodes tess_sizehint=5000, multithreading=false, and uses same tol for initial mesh spacing and final tolerance.

  Better:

  find_candidates(f, region;
      mesh_step = 0.02,
      tolerance = 1e-4,
      tess_sizehint = 50_000,
      maxnodes = 500_000,
      threaded = true,
  )

  Reason: initial mesh spacing controls discovery; tolerance controls candidate location accuracy. Same number for both is often wrong. For
  expensive dispersion_det, multithreading=true may help because RootsAndPoles threads function evaluations.

  2. Add Candidate + Polish Layer
  Do not make solve(region) pretend GRPF roots are final roots. GRPF gives approximate candidate centroids.

  Suggested API:

  cands = find_candidates(plasma, k, region; mesh_step, tolerance)
  roots = polish_candidates(plasma, k, cands; atol=1e-10)

 7. Use Scaled Determinant In Global Search
  GRPF sees f(z) only by phase quadrant, but numerical overflow/underflow still hurts. Use:

  f(ω) = det(D(ω,k)) / prod(row_norms(D))

  or log-safe normalized determinant where possible.

  This makes phase behavior less dominated by huge curl-curl/cold terms and low-frequency scaling.