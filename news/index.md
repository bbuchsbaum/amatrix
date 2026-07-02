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
