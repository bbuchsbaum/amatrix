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

* `matmul`, `crossprod`, `tcrossprod`, `covariance`, `dist_matrix` —
  backend-dispatched with GPU-resident fast paths.
* `chol_factor`, `lu_factor`, `svd_factor` — factorization caching with
  solve, project, and reconstruct methods.
* `eigh`, `am_svd`, `rsvd`, `block_lanczos` — eigendecomposition and
  randomized/iterative spectral methods.
* `sinkhorn` — doubly-stochastic scaling with GPU-resident key-chain loop.
* `ewise`, `am_sweep`, `segment_sum`, `segment_mean` — element-wise and
  reduction operations.
* `woodbury_solve`, `woodbury_logdet` — Woodbury matrix identity solvers.

## Statistical models

* `many_lm` — fit multiple linear models via GPU-accelerated QR with
  coefficient caching.
* `lm_fit`, `array_lm`, `wls_fit` — single and weighted least squares.
* `ridge_path`, `ridge_fit` — ridge regression over lambda grids.
* `lm_loo_cv` — leave-one-out cross-validation via QR downdate.
* `pca_coef` — PCA coefficient extraction from cached SVD factors.

## Sparse operations

* `adgCMatrix` with backend-dispatched SpMV and SpMM via OpenCL, Metal,
  and MLX resident paths.
* `block_lanczos` and `svd_factor` subspace iteration on sparse inputs.

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
