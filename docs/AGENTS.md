# Agent Notes

## Quality Strategy

See `docs/quality-tracking.md` for the full accuracy and performance
methodology. Summary: - **Accuracy:**
[`devtools::test()`](https://devtools.r-lib.org/reference/test.html) —
cross-backend conformance at `1e-4` (GPU) / `1e-10` (CPU) -
**Performance:** `Rscript tools/benchmark-regression.R` — compare to
`tools/baseline.csv` - Every new exported op must be added to the
coverage table in that doc.

## Issue Tracking

Use `bd` (beads) for local issue tracking in this repo.

- The project-local workspace is `.beads/`.
- The issue prefix for this repo is `amatrix`, so new issues must be
  `amatrix-*`.
- After changing issue state, run `bd sync` to commit beads changes to
  git.
- If the workspace looks inconsistent, run `bd doctor` before making
  further changes.
