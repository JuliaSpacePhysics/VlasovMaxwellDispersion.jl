# Development notes

[![Coverage](https://codecov.io/gh/JuliaSpacePhysics/VlasovMaxwellDispersion.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSpacePhysics/VlasovMaxwellDispersion.jl)

## Roadmap

- For `kz = 0`, relativistic paths there need `Im ω ≥ 0`; the cyclotron-maser continuation to damped ω is not implemented.
- Check and potentially reuse `BSplineKit.jl`/`Dierckx.jl`
- Baalrud, S. D. (2013). The incomplete plasma dispersion function: Properties and application to waves in bounded plasmas.
- SciML.HomotopyProblem formalism for continuation

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

