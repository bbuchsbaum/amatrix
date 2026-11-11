# Changelog

## amatrix 0.1.0

Initial public release.

### Core architecture

- S4 classes `adgeMatrix` (dense) and `adgCMatrix` (sparse) extend the
  Matrix package with backend-dispatch slots for transparent GPU
  acceleration.
- Five backend targets: CPU (default), MLX (Apple Silicon), OpenCL
  (CLBlast), Metal (Apple sparse), and ArrayFire (portable GPU).
- Three precision modes: `"strict"` (float64 CPU), `"exact"` (float64
  CPU-pinned), and `"fast"` (float32-permitted GPU).
- Automatic dispatch planning with calibration thresholds, resident
  memory management, and predictable CPU fallback.

### Dense linear algebra

- `matmul`, `crossprod`, `tcrossprod`, `covariance`, `dist_matrix` —
  backend-dispatched with GPU-resident fast paths.
- `chol_factor`, `lu_factor`, `svd_factor` — factorization caching with
  solve, project, and reconstruct methods.
- `eigh`, `am_svd`, `rsvd`, `block_lanczos` — eigendecomposition and
  randomized/iterative spectral methods.
- `sinkhorn` — doubly-stochastic scaling with GPU-resident key-chain
  loop.
- `ewise`, `am_sweep`, `segment_sum`, `segment_mean` — element-wise and
  reduction operations.
- `woodbury_solve`, `woodbury_logdet` — Woodbury matrix identity
  solvers.

### Statistical models

- `many_lm` — fit multiple linear models via GPU-accelerated QR with
  coefficient caching.
- `lm_fit`, `array_lm`, `wls_fit` — single and weighted least squares.
- `ridge_path`, `ridge_fit` — ridge regression over lambda grids.
- `lm_loo_cv` — leave-one-out cross-validation via QR downdate.
- `pca_coef` — PCA coefficient extraction from cached SVD factors.

### Sparse operations

- `adgCMatrix` with backend-dispatched SpMV and SpMM via OpenCL, Metal,
  and MLX resident paths.
- `block_lanczos` and `svd_factor` subspace iteration on sparse inputs.

### Backend packages (separate installation)

- `amatrix.mlx` — Apple MLX backend via mlx-c; full dense + sparse
  support.
- `amatrix.opencl` — OpenCL + CLBlast backend; dense BLAS,
  factorizations, and sparse products.
- `amatrix.metal` — Apple Metal backend; sparse-first with MPSGraph
  kernels.
- `amatrix.arrayfire` — ArrayFire backend; portable GPU with safety
  gates for Apple Silicon.

### Benchmark harness

- `tools/benchmark-regression.R` — canonical regression harness with
  per-backend worker isolation, direct-file and source-entry parity, and
  automatic baseline comparison.

### Observability and honest defaults

- [`amatrix_backend_health_probe()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_health_probe.md)
  — run a canary op against a backend and record its health (`healthy` /
  `unhealthy:<reason>`). Every backend is probed on first use;
  subsequent routing decisions respect the probe.
- [`amatrix_backend_status()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_status.md)
  now reports per-backend `health` and `health_reason` columns alongside
  capabilities, features, and precision modes.
- [`amatrix_fallback_log()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_fallback_log.md)
  and
  [`amatrix_fallback_log_reset()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_fallback_log_reset.md)
  — structured log of every dispatch fall-through from a backend to CPU.
  A non-empty log after a clean conformance run is a stop-ship
  condition.
- [`amatrix_benchmark_report()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_benchmark_report.md)
  — reads `tools/baseline.csv` and the cached calibration, returning
  cold-vs-warm timings, warm/cold ratios, calibrated thresholds, and
  `gpu_wins` flags per (op × backend).
- Calibration cache is tagged with a `sys_hash` over OS, machine,
  platform, and R version. Stale caches from a different machine are
  rejected automatically instead of producing wrong routing.
- Auto-fallback with telemetry: every backend error in
  `amatrix_dispatch_op` is re-signaled as a calling-style condition (so
  `withCallingHandlers` can observe the original error class) AND emits
  an `amatrix_fallback` condition with structured metadata, before
  routing to CPU.
- `mode = "balanced"` is now deprecated. It was never fully implemented
  (routed to CPU under the hood) and now maps to `"exact"` with a
  one-time-per-session deprecation warning. Use `"exact"` or `"fast"`.

### Pre-release review fixes

- Drop-in semantics: `x[i]` / `x[i] <- v` index elements as in base R;
  `log(x, base)` honours `base`;
  [`norm()`](https://rdrr.io/r/base/norm.html) matches
  [`base::norm()`](https://rdrr.io/r/base/norm.html) for every type and
  no longer alters [`norm()`](https://rdrr.io/r/base/norm.html) for
  plain matrices; [`t()`](https://rdrr.io/r/base/t.html) views support
  the usual matrix operations; results keep their amatrix class across
  sub-assignment, products with vectors and plain matrices, and
  [`kronecker()`](https://rdrr.io/r/base/kronecker.html).
- Identity-keyed caches can no longer return a factor of different data
  (after `-x`, [`round()`](https://rdrr.io/r/base/Round.html),
  [`signif()`](https://rdrr.io/r/base/Round.html),
  [`zapsmall()`](https://rdrr.io/r/base/zapsmall.html), or re-used
  weights), and cache entries are released when their source matrix is
  collected.
- Loading the package, probing a backend and `trace_estim(seed = )` no
  longer touch the user’s random-number stream.
- [`block_lanczos()`](https://bbuchsbaum.github.io/amatrix/reference/block_lanczos.md)
  /
  [`block_svd()`](https://bbuchsbaum.github.io/amatrix/reference/block_lanczos.md)
  check residuals against a new `tol` argument, extend the iteration
  automatically and warn (`amatrix_not_converged`) when the tolerance is
  not met.
- [`lm_fit()`](https://bbuchsbaum.github.io/amatrix/reference/lm_fit.md)
  /
  [`wls_fit()`](https://bbuchsbaum.github.io/amatrix/reference/wls_fit.md)
  raise `amatrix_rank_deficient` for a singular design and warn
  (`amatrix_ill_conditioned`) when the normal equations are unreliable;
  [`lm_loo_cv()`](https://bbuchsbaum.github.io/amatrix/reference/lm_loo_cv.md)
  uses the exact leverage identity;
  [`lu_factor()`](https://bbuchsbaum.github.io/amatrix/reference/lu_factor.md)
  stores a real LU factorization.
- [`amatrix_calibrate()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_calibrate.md)
  writes to the user cache directory only with `persist = TRUE`.
- `dist_matrix(method = "cosine")` returns the cosine distance
  (`1 - similarity`); the similarity is
  `kernel_matrix(kernel = "cosine")`.
- [`sweep()`](https://rdrr.io/r/base/sweep.html) now dispatches on
  amatrix objects (it is exported as a standard generic with
  [`base::sweep`](https://rdrr.io/r/base/sweep.html) as default);
  [`am_sweep()`](https://bbuchsbaum.github.io/amatrix/reference/am_sweep.md)
  defaults to `FUN = "-"` like base;
  [`rowscale()`](https://bbuchsbaum.github.io/amatrix/reference/rowscale.md),
  [`colscale()`](https://bbuchsbaum.github.io/amatrix/reference/rowscale.md)
  and multiplicative sweeps keep sparse input sparse.
  [`mean()`](https://rdrr.io/r/base/mean.html),
  [`which()`](https://rdrr.io/r/base/which.html) and
  [`isSymmetric()`](https://rdrr.io/r/base/isSymmetric.html) work with
  only amatrix attached. `dim(x) <- value` returns a valid, freshly
  identified object.
- Much lower per-operation overhead and fewer data copies on the CPU
  path.
- Attaching the package no longer masks other packages’ functions: the
  matrix trace is
  [`mat_trace()`](https://bbuchsbaum.github.io/amatrix/reference/mat_trace.md)
  (was `trace`, which masked
  [`base::trace`](https://rdrr.io/r/base/trace.html)), the truncated-SVD
  wrapper is
  [`am_irlba()`](https://bbuchsbaum.github.io/amatrix/reference/am_irlba.md)
  (was `irlba`, which masked `irlba::irlba`), and `svd` is a standard
  generic with [`base::svd`](https://rdrr.io/r/base/svd.html) as its
  default.

### Error handling and diagnostics

- User-facing errors in wrapper functions now emit classed conditions
  (`amatrix_bad_arg`, `amatrix_bad_backend`, `amatrix_backend_exists`,
  `amatrix_subspace_error`) with `call = NULL` so messages don’t leak
  internal call stacks. `tryCatch(..., amatrix_bad_arg = ...)` and
  `expect_error(class = "amatrix_bad_arg")` both work.

### Quality machinery

- `planning_docs/quality-tracking.md` is the authoritative quality doc.
  It carries the coverage matrix (every exported op × four test types +
  benchmark row + per-backend tier), the stop-ship rules, and the honest
  backend tier assessment.
- `tools/audit-dispatch.R` — static AST-based audit that enumerates
  every required dispatch signature for `%*%` / `crossprod` /
  `tcrossprod` on mixed plain/amatrix pairs and fails the PR gate if any
  is missing.
- Residency tripwire (`options(amatrix.residency.tripwire = TRUE)` or
  `AMATRIX_RESIDENCY_TRIPWIRE=1`) increments a counter on every real
  GPU-to-host copy; the conformance suite asserts zero trips after a
  clean run.
- Three-lane CI: PR gate (`.github/workflows/R-CMD-check.yaml` +
  fast-gates job), nightly stress (`nightly-stress.yaml`), release gate
  (`release-gate.yaml` with required reviewers environment).
