# amatrix 0.1.0

Initial public release.

## Core architecture

* S4 classes `adgeMatrix` (dense) and `adgCMatrix` (sparse) extend the
  Matrix package with backend-dispatch slots for transparent GPU acceleration.
* Five backend targets: CPU (default), MLX (Apple Silicon), OpenCL (CLBlast),
  Metal (Apple sparse), and ArrayFire (portable GPU).
* Three precision modes: `"strict"` (float64 CPU), `"exact"` (float64
  CPU-pinned), and `"fast"` (float32-permitted GPU).
* Automatic dispatch planning with calibration thresholds, resident memory
  management, and predictable CPU fallback.

## Dense linear algebra

* `am_matmul`, `crossprod`, `tcrossprod`, `covariance`, `am_dist_matrix` —
  backend-dispatched with GPU-resident fast paths.
* `am_chol_factor`, `am_lu_factor`, `am_svd_factor` — factorization caching with
  solve, project, and reconstruct methods.
* `am_eigh`, `am_svd`, `rsvd`, `am_block_lanczos` — eigendecomposition and
  randomized/iterative spectral methods.
* `am_sinkhorn` — doubly-stochastic scaling with GPU-resident key-chain loop.
* `am_ewise`, `am_sweep`, `am_segment_sum`, `am_segment_mean` — element-wise and
  reduction operations.
* `am_woodbury_solve`, `am_woodbury_logdet` — Woodbury matrix identity solvers.

## Statistical models

* `am_many_lm` — fit multiple linear models via GPU-accelerated QR with
  coefficient caching.
* `am_lm_fit`, `am_array_lm`, `am_wls_fit` — single and weighted least squares.
* `am_ridge_path`, `am_ridge_fit` — ridge regression over lambda grids.
* `am_lm_loo_cv` — leave-one-out cross-validation via QR downdate.
* `am_pca_coef` — PCA coefficient extraction from cached SVD factors.

## Sparse operations

* `adgCMatrix` with backend-dispatched SpMV and SpMM via OpenCL, Metal,
  and MLX resident paths.
* `am_block_lanczos` and `am_svd_factor` subspace iteration on sparse inputs.

## Backend packages (separate installation)

* `amatrix.mlx` — Apple MLX backend via mlx-c; full dense + sparse support.
* `amatrix.opencl` — OpenCL + CLBlast backend; dense BLAS, factorizations,
  and sparse products.
* `amatrix.metal` — Apple Metal backend; sparse-first with MPSGraph kernels.
* `amatrix.arrayfire` — ArrayFire backend; portable GPU with safety gates
  for Apple Silicon.

## Benchmark harness

* `tools/benchmark-regression.R` — canonical regression harness with
  per-backend worker isolation, direct-file and source-entry parity, and
  automatic baseline comparison.

## Observability and honest defaults

* `amatrix_backend_health_probe()` — run a canary op against a backend and
  record its health (`healthy` / `unhealthy:<reason>`). Every backend is
  probed on first use; subsequent routing decisions respect the probe.
* `amatrix_backend_status()` now reports per-backend `health` and
  `health_reason` columns alongside capabilities, features, and precision
  modes.
* `amatrix_fallback_log()` and `amatrix_fallback_log_reset()` — structured
  log of every dispatch fall-through from a backend to CPU. A non-empty log
  after a clean conformance run is a stop-ship condition.
* `amatrix_benchmark_report()` — reads `tools/baseline.csv` and the cached
  calibration, returning cold-vs-warm timings, warm/cold ratios, calibrated
  thresholds, and `gpu_wins` flags per (op × backend).
* Calibration cache is tagged with a `sys_hash` over OS, machine, platform,
  and R version. Stale caches from a different machine are rejected
  automatically instead of producing wrong routing.
* Auto-fallback with telemetry: every backend error in `amatrix_dispatch_op`
  is re-signaled as a calling-style condition (so `withCallingHandlers` can
  observe the original error class) AND emits an `amatrix_fallback`
  condition with structured metadata, before routing to CPU.
* `mode = "balanced"` is now deprecated. It was never fully implemented
  (routed to CPU under the hood) and now maps to `"exact"` with a
  one-time-per-session deprecation warning. Use `"exact"` or `"fast"`.

## Pre-release review fixes

* Drop-in semantics: `x[i]` / `x[i] <- v` index elements as in base R;
  `log(x, base)` honours `base`; `norm()` matches `base::norm()` for every
  type and no longer alters `norm()` for plain matrices; `t()` views support
  the usual matrix operations; results keep their amatrix class across
  sub-assignment, products with vectors and plain matrices, and `kronecker()`.
* Identity-keyed caches can no longer return a factor of different data
  (after `-x`, `round()`, `signif()`, `zapsmall()`, or re-used weights), and
  cache entries are released when their source matrix is collected.
* Loading the package, probing a backend and `am_trace_estim(seed = )` no longer
  touch the user's random-number stream.
* `am_block_lanczos()` / `am_block_svd()` check residuals against a new `tol`
  argument, extend the iteration automatically and warn
  (`amatrix_not_converged`) when the tolerance is not met.
* `am_lm_fit()` / `am_wls_fit()` raise `amatrix_rank_deficient` for a singular
  design and warn (`amatrix_ill_conditioned`) when the normal equations are
  unreliable; `am_lm_loo_cv()` uses the exact leverage identity; `am_lu_factor()`
  stores a real LU factorization.
* `amatrix_calibrate()` writes to the user cache directory only with
  `persist = TRUE`.
* `am_dist_matrix(method = "cosine")` returns the cosine distance
  (`1 - similarity`); the similarity is `am_kernel_matrix(kernel = "cosine")`.
* `sweep()` now dispatches on amatrix objects (it is exported as a standard
  generic with `base::sweep` as default); `am_sweep()` defaults to `FUN = "-"`
  like base; `am_rowscale()`, `am_colscale()` and multiplicative sweeps keep sparse
  input sparse. `mean()`, `which()` and `isSymmetric()` work with only
  amatrix attached. `dim(x) <- value` returns a valid, freshly identified
  object.
* Much lower per-operation overhead and fewer data copies on the CPU path.
* Accelerator path: with `options(amatrix.defer_host = TRUE)` a successful
  dense product, cross-product or element-wise result was mistaken for a
  failure, which then marked the backend unhealthy and silently moved all
  work to the CPU; device-resident chains now work. When uncalibrated, small
  matrices and memory-bound operations on host data stay on the CPU instead
  of paying for a device round trip. The planner no longer resolves the
  calibration cache path on every operation.
* For package authors: a new article, `vignette("package-authors")`; a
  testing kit (`amatrix_test_backend()`, `amatrix_with_test_backend()`,
  `amatrix_skip_if_no_accelerator()`, `amatrix_has_accelerator()`); and
  `mode = "float32"` as an explicit spelling of `mode = "fast"`.
* Naming: the exported API uses two prefixes. `amatrix_*` covers session,
  backend and diagnostic functions; `am_*` covers everything that computes
  (`am_matmul()`, `am_chol_factor()`, `am_svd_factor()`, `am_many_lm()`,
  `am_lm_fit()`, `am_rsvd()`, `am_irlba()`, `am_trace()`, `am_sqrtm()`, ...).
  Class constructors and coercions (`adgeMatrix()`, `adgCMatrix()`,
  `as_adgeMatrix()`, `kron_matrix()`, `resident_handle()`) and
  `with_amatrix()` keep their names. Attaching the package masks nothing from
  base R, the recommended packages, or packages such as 'irlba', 'rsvd' and
  'rlang'; `svd`, `norm`, `sweep`, `mean`, `which` and `isSymmetric` are
  exported only as standard generics whose default is the base function.

## Error handling and diagnostics

* User-facing errors in wrapper functions now emit classed conditions
  (`amatrix_bad_arg`, `amatrix_bad_backend`, `amatrix_backend_exists`,
  `amatrix_subspace_error`) with `call = NULL` so messages don't leak
  internal call stacks. `tryCatch(..., amatrix_bad_arg = ...)` and
  `expect_error(class = "amatrix_bad_arg")` both work.

## Quality machinery

* `planning_docs/quality-tracking.md` is the authoritative quality doc. It
  carries the coverage matrix (every exported op × four test types +
  benchmark row + per-backend tier), the stop-ship rules, and the honest
  backend tier assessment.
* `tools/audit-dispatch.R` — static AST-based audit that enumerates every
  required dispatch signature for `%*%` / `crossprod` / `tcrossprod` on
  mixed plain/amatrix pairs and fails the PR gate if any is missing.
* Residency tripwire (`options(amatrix.residency.tripwire = TRUE)` or
  `AMATRIX_RESIDENCY_TRIPWIRE=1`) increments a counter on every real
  GPU-to-host copy; the conformance suite asserts zero trips after a
  clean run.
* Three-lane CI: PR gate (`.github/workflows/R-CMD-check.yaml` +
  fast-gates job), nightly stress (`nightly-stress.yaml`), release gate
  (`release-gate.yaml` with required reviewers environment).
