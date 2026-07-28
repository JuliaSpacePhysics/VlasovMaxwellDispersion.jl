# CLAUDE.md

Keep in-repo memory (project knowledge) at `.claude/memory/`: `MEMORY.md` as index.

## Development status

Under active development, not yet stable.

- **Breaking changes allowed and encouraged** when they improve the design; no
  backward compatibility for its own sake.
- **Prefer the cleanest API**. Refactor freely; update tests and docs to match.

## Layout & workflow

Tests are `@testitem`s under `test/{distributions,solvers,kernels,relativistic}/`

- Focused: `julia --project=test -e 'using TestItemRunner;
  TestItemRunner.run_tests(pwd(); filter = ti -> occursin(r"...", ti.name))'`.
- `test/standalone/` deliberately kept out of the suite
- Full runs are slow; reach only when the blast radius is genuinely package-wide.
