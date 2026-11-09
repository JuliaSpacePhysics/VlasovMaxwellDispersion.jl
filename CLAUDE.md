# CLAUDE.md

## Development status

Under active development, not yet stable.

- **Breaking changes allowed and encouraged** when they improve the design.
- **Prefer clean API**. Refactor freely; update tests and docs to match.

## Layout & workflow

- `test/standalone/` deliberately kept out of the suite
- Full runs are slow; reach only when the blast radius is genuinely package-wide.

## Performance work

- A solve costs (χ evaluations) × (cost per evaluation), and a plan change can move them in opposite directions. Judge by the global solve benchmarks and a case-study survey (`nevals`, retcode, root set), not by `plan`/`evaluate` alone.
- Table accuracy (`FixedNodeEval.tol`) is the lever, not solver tolerances: loosening AAA to a table's error may drop real weak roots while reporting `Success`.
- Compare survey root sets by matching (k, ω) points, not branch counts;
  `proton_shell_Min15` is `Saturated` and its set varies between plan variants.
