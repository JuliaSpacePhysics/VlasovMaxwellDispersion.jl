# CLAUDE.md

## Agent memory

In-repo persistent memory (project knowledge) is at `.claude/memory/`: one fact per file with frontmatter, `MEMORY.md` as index. 
Read index before starting work; save/update/link memories there.

## Development status

Under active development, not yet stable.

- **Breaking changes allowed and encouraged** when they improve the design; no
  backward compatibility for its own sake.
- **Prefer the cleanest API**. Refactor freely; update tests and docs to match.

## Layout & workflow

- Tests are `@testitem`s: `julia --project=test -e 'using TestItemRunner;
  TestItemRunner.run_tests(pwd(); filter = ti -> occursin(r"...", ti.name))'`.
