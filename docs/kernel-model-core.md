# Kernel And Model-Core Notes

This note is for package authors and for future `amatrix` extension packages.

The current architecture has three layers:

- object layer
  - `adgeMatrix`, `adgCMatrix`
- kernel layer
  - `am_matmul()`, `am_crossprod()`, `am_tcrossprod()`, `am_solve()`, `am_ewise()`, `am_row_sums()`, `am_col_sums()`
- model-core layer
  - `am_lm_fit()`
  - `am_ridge_fit()`

## What To Call

Use:

- public generics for the user-facing API surface
  - examples: `%*%`, `crossprod(x)`, `qr(x)` when that path exists
- internal `am_*` wrappers for package implementation paths
  - examples: `am_matmul(x, y)`, `am_crossprod(x, y)`, `am_solve(x, b)`

That split matters because `amatrix` internals should not depend on ambient dispatch being cooperative.

## Current Stable Kernel Surface

These wrappers are the current intended implementation surface:

- `am_matmul(x, y)`
- `am_crossprod(x, y = NULL)`
- `am_tcrossprod(x, y = NULL)`
- `am_solve(a, b = NULL)`
- `am_ewise(op, e1, e2 = NULL)`
- `am_row_sums(x, na.rm = FALSE, dims = 1L)`
- `am_col_sums(x, na.rm = FALSE, dims = 1L)`

The wrappers:

- preserve `amatrix` object semantics where appropriate
- go through backend planning instead of backend-name conditionals
- respect `strict` vs `fast` precision policy

## Current Model-Core Surface

The first model-core entry points are:

- `am_lm_fit(X, Y, ...)`
- `am_ridge_fit(X, Y, lambda, ...)`
- `am_many_lm(X, Y, ...)`
- `am_array_lm(X, Y, ...)`
- `am_covariance(X, ...)`
- `am_correlation(X, ...)`

These functions are designed around the highest-value workload in the PRD:

- one shared design matrix `X`
- many response matrices or vectors `Y`

Current center of gravity:

- `am_many_lm(X, Y, method = "qr")`

That is the strongest current surface for package authors who want to swap in a small number of lines and get a real workload-level speedup.

They currently provide:

- Matrix-compatible input handling
- backend-neutral kernel use
- internal shared-`X` cache reuse for repeated fits
- a repeated-response workflow surface that exposes response counts and per-response residual summaries
- an array-aware regression surface that can restore fitted and residual outputs to response-array shape
- covariance and correlation helpers that route the heavy second-moment step through `am_crossprod()`

## Shared-X Cache Rule

The cache is intentionally narrow:

- it is internal-only
- it is keyed by stable `amatrix` identity plus backend/policy/precision
- it currently reuses `crossprod(X)` and rank metadata
- it is never required for correctness

This is not a general user-visible state mechanism. It is a workload-specific optimization.

## What Not To Do

Do not:

- call backend-native bridges directly from higher-level code
- scatter `if (backend == "mlx")` logic through algorithms
- assume resident execution is available
- assume `fast` precision unless the caller opts in

## Package Structure Guidance

Do not create one package per upstream CRAN algorithm clone by default.

Prefer:

- `amatrix`
  - substrate and kernels
- `amatrix.models`
  - fitters and repeated-solve workloads
- `amatrix.algorithms`
  - randomized SVD, truncated methods, PCA helpers

That keeps the project looking like a numerical platform rather than a shadow collection of GPU-flavored package forks.
