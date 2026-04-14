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

This layering is deliberate. `amatrix` should look like infrastructure, not a catalog of standalone GPU wrappers for individual statistics functions.

## What To Call

Use:

- public generics for the user-facing API surface
  - examples: `%*%`, `crossprod(x)`, `qr(x)` when that path exists
- internal `am_*` wrappers for package implementation paths
  - examples: `am_matmul(x, y)`, `am_crossprod(x, y)`, `am_solve(x, b)`

That split matters because `amatrix` internals should not depend on ambient dispatch being cooperative.

## Product Syntax Boundary

Use the operator and kernel surfaces deliberately:

- ordinary products
  - `%*%`
- transpose-heavy products
  - `crossprod(x, y)`
  - `tcrossprod(x, y)`
- full DGEMM control
  - internal `.am_gemm(...)` today
  - planned public `am_gemm(A, B, C = NULL, alpha = 1, beta = 0, transA = FALSE, transB = FALSE)`

Today `crossprod()` and `tcrossprod()` are the reliable hot-path forms for
transpose algebra because `t(adgeMatrix)` still uses a materialized stepping
stone rather than a dedicated structural view. The current `src_id` shortcut
fixes the GPU re-upload problem, but it still pays host transpose work. A
Matrix-compatible transpose view is the intended next step for `t(A) %*% B` and
`A %*% t(B)`.

What is not in scope:

- no general lazy expression tree behind `%*%`, `+`, and `*`
- no requirement that eager operator chains such as `A %*% B + C` fuse
  automatically

## Current Stable Kernel Surface

These wrappers are the current intended implementation surface:

**Matrix products and reductions:**
- `am_matmul(x, y)`
- `am_crossprod(x, y = NULL)`
- `am_tcrossprod(x, y = NULL)`
- `am_solve(a, b = NULL)`
- `am_ewise(op, e1, e2 = NULL)`
- `am_row_sums(x, na.rm = FALSE, dims = 1L)`
- `am_col_sums(x, na.rm = FALSE, dims = 1L)`

**Planned explicit product-control kernel:**
- `am_gemm(A, B, C = NULL, alpha = 1, beta = 0, transA = FALSE, transB = FALSE)`
  - explicit BLAS-style product control for package authors who need fused
    accumulator semantics

**Factorizations and decompositions (GPU-accelerated where available):**
- `am_qr(x)` / `am_qr_factor(x)` → `amQR`
- `am_chol(x)` / `am_chol_factor(x)` → `amChol`
- `am_chol_solve(factor, B)` — batched forward/backward substitution
- `am_rsvd(x, k, n_iter, n_extra)` — fast-mode backend randomized SVD with CPU fallback under strict precision
- `am_svd_factor(x, k)` → `amSVD` — cached SVD factor
- `am_svd_project(factor, Y)` — batched projection into latent space
- `am_svd_reconstruct(factor, scores)` — back-projection
- `am_block_lanczos(x, nv, ...)` — GEMM-oriented truncated SVD prototype for dense MLX workloads
- `am_irlba(x, k, ...)` — compatibility Lanczos wrapper (GPU matvecs; use `am_rsvd` or `am_block_lanczos` for speed)
- `am_irlba_native(x, k, ...)` — ArrayFire-native Lanczos
- `am_eigen(x)` — eigendecomposition (CPU today; GPU syevd is planned)

**Similarity and distance:**
- `am_dist(x, y, method)` — pairwise distance matrix (GPU via tcrossprod)
- `am_kernel(x, y, kernel, ...)` — pairwise kernel matrix (GPU bridge)
- `am_covariance(x, ...)` — sample covariance
- `am_correlation(x, ...)` — sample correlation

The wrappers:

- preserve `amatrix` object semantics where appropriate
- go through backend planning instead of backend-name conditionals
- respect `strict` vs `fast` precision policy

## Specialized Intermediates

Specialized S4 intermediates are allowed when they preserve matrix semantics
and avoid unnecessary materialization.

Current examples:

- `amQR`
- `amChol`
- `amSVD`
- the current `src_id`-linked transposed `adgeMatrix` shortcut

Planned structural example:

- `aTransposeView` or a broader `aMatrixView` family that replaces the current
  transpose shortcut with a true structural view

These are sanctioned performance surfaces. They are not a general deferred
evaluation DSL.

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
- it currently reuses `crossprod(X)`, QR factors, and rank metadata
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
  - substrate, kernels, backend registry, dispatch
  - all `am_*` wrappers including decompositions and distance/kernel
- `amatrix.models`
  - repeated-response least squares (`many_lm`, `array_lm`)
  - ridge, weighted least squares
  - covariance and correlation helpers
- `amatrix.algorithms`
  - `rsvd()` — randomized SVD (wraps `am_rsvd`, adds cross-validation, sketch reuse)
  - `pca()` — GPU-accelerated PCA via `am_svd_factor` + projection helpers
  - `block_lanczos()` — GPU-native truncated SVD via block Lanczos
  - `eigensolver()` — GPU symmetric eigendecomp when implemented
  - `pcr()` — principal components regression composing `am_svd_factor` + OLS
  - `cca()` — canonical correlation via paired SVD factors
- `amatrix.image` (future)
  - neuroimaging and image-domain kernels

That keeps the project looking like a numerical platform rather than a shadow collection of GPU-flavored package forks.

## GPU Acceleration Status by Kernel

| Kernel | CPU | MLX GPU | ArrayFire GPU | Next step |
|--------|-----|---------|---------------|-----------|
| matmul / crossprod / tcrossprod | ✅ | ✅ fast | ✅ fast | — |
| QR | ✅ | ✅ fast (resident + compact TSQR) | ✅ fast | resident-Q is the premium QR path; compact TSQR is acceptable when factor shape matters |
| rsvd | ✅ irlba / svdr fallback | ✅ native fast path | ⚠️ wrapper path present | decide ArrayFire product status and keep the fast-path contract explicit |
| block Lanczos / truncated SVD | ✅ irlba / svdr | ⚠️ prototype fast path | ❌ | decide promotion criteria; current MLX prototype now beats CPU at `3000x1200` and `5000x2000`, with mildly oversampled defaults and an explicit ~4-5% approximation envelope |
| Cholesky | ✅ | ✅ bridge-backed factor path | ⚠️ CPU via stub | keep ArrayFire fallback boundary explicit; broaden adoption guidance as needed |
| trsm (batched solve) | ✅ | ✅ bridge-backed triangular solve | ⚠️ sequential / CPU solve | validated on representative SPD workflows; next frontier is elsewhere |
| Eigendecomp (symmetric) | ✅ base::eigen | ❌ | ❌ | GPU syevd |
| Distance / kernel | ✅ | ✅ tcrossprod | ✅ AF bridge | tiled large-n |
| Sparse × dense | ✅ Matrix | ❌ | ❌ | deferred until after product hardening |
