# Flagship Workflows

This note shows the current center of gravity for `amatrix`.

The project is not trying to be "GPU for every matrix call." The strongest current public story is:

- one shared design matrix `X`
- many response columns `Y`
- cached repeated-fit work
- QR-backed least squares
- backend acceleration when explicitly allowed and numerically appropriate

The current flagship entry point is:

```r
library(amatrix.models)

fit <- many_lm(X, Y_many, method = "qr", cache = TRUE)
```

This is where the current QR work is paying off most clearly.

## Algebra Boundary

The flagship benchmark story is not "turn every R operator chain into a lazy
tensor graph." The intended boundary is:

- ordinary products
  - `%*%`
- transpose-heavy products
  - `crossprod()` and `tcrossprod()` today
  - a proper transpose view next, replacing the current `src_id` shortcut
- full BLAS-style control
  - planned explicit `am_gemm()`

That keeps the surface honest while leaving the real speed fight focused on
cached QR, Cholesky, and shared-`X` many-`Y` workloads.

## Current Flagship

Use `amatrix.models::many_lm()` with `method = "qr"` when the workload is one `X` and many response columns:

```r
X <- adgeMatrix(design, mode = "fast", backend = "mlx")
fit <- many_lm(X, Y_many, method = "qr", cache = TRUE)

fit$responses
fit$rss
fit$sigma2
coef(fit)
```

Why this is the flagship:

- the shared-`X` cache is reused automatically
- the QR path is numerically stronger than normal equations
- the resident MLX QR path now starts to dominate once the number of right-hand sides grows

Current many-RHS benchmark note on this machine for `1024x128`:

| rhs_cols | base cached QR | MLX native resident | speedup |
|---|---|---|---|
| 8   | 0.0016 s | 0.0020 s | 0.8x (CPU wins) |
| 32  | 0.0058 s | 0.0018 s | 3.2x |
| 128 | 0.0215 s | 0.0025 s | **8.6x** |

At 128 right-hand sides the MLX resident path is **8.6× faster** than the CPU cached-QR baseline. This is the headline number for the shared-X many-Y use case.

## Second Flagship: GPU Cholesky + Batched Solve

The Cholesky pattern is the second highest-impact surface: factorize a symmetric positive-definite matrix once, then solve for many right-hand sides.

```r
library(amatrix)

K <- adgeMatrix(
  am_kernel(X, X, kernel = "rbf", sigma = 1.0),
  mode = "fast",
  backend = "mlx"
)
L <- am_chol_factor(K)   # GPU-backed Cholesky factor

# Solve K %*% alpha = y for many y columns — one batched trsm call
alpha <- am_chol_solve(L, Y_many)
```

Why this matters:
- Ridge regression: factorize `(X^TX + λI)` once, solve for many λ values or Y columns
- Gaussian process: Cholesky of kernel matrix, predict at many test points
- lme4-style mixed models: Cholesky of precision matrix, many PIRLS solve iterations
- Mahalanobis distance: `(x - μ)^T Σ^{-1} (x - μ)` via Cholesky of Σ

Current benchmark note on this machine, using `Rscript -e 'source("tools/benchmark-cholesky-runtime.R", local = TRUE)'`:

- ridge-like SPD, `768x768`, `rhs_cols = 64`
  - CPU factor: about `0.065 s`
  - MLX factor: about `0.019 s`
  - CPU factor + batched solve: about `0.080 s`
  - MLX factor + batched solve: about `0.029 s`
- kernel-like SPD, `640x640`, `rhs_cols = 32`
  - CPU factor: about `0.036 s`
  - MLX factor: about `0.011 s`
  - CPU factor + batched solve: about `0.042 s`
  - MLX factor + batched solve: about `0.017 s`

Quality stayed tight in the same run:

- factor reconstruction residuals: about `2e-7`
- solve error versus CPU reference: about `5e-7`

Status: validated on native MLX for representative ridge-like and kernel-like SPD many-RHS workloads. ArrayFire is still CPU fallback/stub territory for Cholesky and triangular solve on this workflow.

---

## Third Flagship: GPU-Native PCA via am_rsvd

For large matrices where full SVD is infeasible, `am_rsvd` provides a GPU-native randomized SVD with a single device-side evaluation:

```r
library(amatrix)

X <- adgeMatrix(data_matrix, mode = "fast", backend = "mlx")

# GPU-native randomized SVD — zero per-step CPU syncs
svd_result <- am_rsvd(X, k = 50, n_iter = 3)

# Cache the factor for repeated projections
fac <- am_svd_factor(X, k = 50)
scores  <- am_svd_project(fac, Y_new)     # project new data
recon   <- am_svd_reconstruct(fac, scores) # reconstruct
```

Why this matters:
- fMRI: SVD of (time × voxels) matrix, project many condition contrasts
- GWAS: LD-adjusted regression via SVD of genotype matrix
- Regularization path: vary λ on fixed SVD without refactoring
- Any PCA/reduced-rank workflow replacing `prcomp()` or `irlba::irlba()`

Performance: beats LAPACK rsvd at m,n≥3000 with k≤100. Does NOT beat CPU for small matrices — dispatch threshold enforced.

Status: **Implemented** in MLX. ArrayFire exposes an `rsvd` path as well, but quality/performance validation is still a separate hardening task.

---

## Three Surfaces

### Many Responses
This is the primary surface today. It is the one to try first if your package or workflow has:

- one design matrix
- many response vectors or columns
- repeated fit/reuse structure

### Array Responses

Use `amatrix.models::array_lm()` when the responses are naturally array-shaped and you want fitted values or residuals restored to that shape:

```r
X <- adgeMatrix(design)
fit <- array_lm(
  X,
  Y_array,
  weights = w,
  method = "qr",
  include_fitted = TRUE,
  include_residuals = TRUE
)

dim(fitted(fit))
dim(residuals(fit))
fit$response_dims
```

This is still a general array-response API. It is not branded around any one scientific domain.

### Single Fit

Use `amatrix.models::lm_fit()` when you want one least-squares fit object with explicit method choice:

```r
library(amatrix)
library(amatrix.models)

X <- adgeMatrix(design)
fit <- lm_fit(X, Y, method = "qr")
coef(fit)
fitted(fit)
residuals(fit)
```

Current methods:

- `"normal"`
  - normal equations via `crossprod(X)` and `solve()`
- `"qr"`
  - QR-backed least squares via `am_qr()` and QR helper methods

### Weighted Fits

Use `amatrix.models::wls_fit()` when you want a single weighted fit object:

```r
fit <- wls_fit(X, Y, weights = w, method = "qr")
coef(fit)
```

The current implementation uses row-weighted transformed design and response matrices, then routes through the same shared model-core machinery. The same weighting path now also underlies `many_lm(..., weights = w)` and `array_lm(..., weights = w)`.

### Covariance And Correlation

Use `amatrix.models::covariance()` and `amatrix.models::correlation()` when the main task is matrix-to-matrix similarity structure rather than repeated regression:

```r
S <- covariance(X)
R <- correlation(X)
S_block <- covariance(X, block_size = 256L)
```

These helpers currently:

- accept ordinary matrix-like or `amatrix` inputs
- keep the heavy second-moment step on `am_crossprod()`
- return ordinary `adgeMatrix` results
- support weighted covariance through `weights = w`
- support blockwise evaluation through `block_size` when you want to limit the width of each second-moment multiply

Current benchmark note on this machine:

- `covariance`: about `0.0143 s`
- `weighted_covariance`: about `0.0170 s`
- `correlation`: about `0.0150 s`

## Mode And Backend

The constructor `adgeMatrix(x, mode=, backend=)` is the primary user API.

| mode | semantics |
|---|---|
| `"exact"` | strict float64, CPU-pinned — no GPU, no silent downcast |
| `"balanced"` | strict float64, currently conservative CPU semantics unless an explicit backend is supplied |
| `"fast"` | float32-permitted execution; specify `backend =` to activate an accelerator path |

Default (`mode` omitted) uses `"balanced"` semantics: strict precision, CPU unless a backend is specified.

For accelerator-oriented execution, specify both `mode` and `backend`:

```r
X <- adgeMatrix(design, mode = "fast", backend = "mlx")
fit <- am_many_lm(X, Y_many, method = "qr", cache = TRUE)
```

The `backend=` argument is currently the explicit accelerator opt-in. `balanced` is kept conservative until automatic routing is hardened enough to document confidently.

## Why Cache Matters

The current shared-`X` cache is internal and narrow:

- keyed by stable `amatrix` identity
- reuses `crossprod(X)` for normal-equation fits
- reuses `am_qr(X)` for QR-backed fits
- never changes object semantics

On this machine, the current benchmark note shows:

- `am_many_lm(..., method = "normal")`
  - `cache_off`: about `0.0733 s`
  - `cache_on`: about `0.0390 s`
- `am_many_lm(..., method = "qr")`
  - `cache_off`: about `0.0637 s`
  - `cache_on`: about `0.0287 s`

See [backend-benchmarks.md](/Users/bbuchsbaum/code/amatrix/docs/backend-benchmarks.md) for the current numbers.

## Current Recommendation

If you are choosing among the current public surfaces:

- use `amatrix.models::many_lm(..., method = "qr")` first when the workload is one `X`, many response columns
- use `amatrix.models::array_lm()` when the response structure is naturally array-shaped, optionally weighted
- use `amatrix.models::lm_fit()` for one explicit fit
- use `amatrix.models::wls_fit()` when you want a single weighted fit object
- use `amatrix.models::covariance()` or `amatrix.models::correlation()` for similarity-structure workloads

That is the most mature part of the project today.

---

## Reproducing These Benchmarks

All benchmark scripts live in `tools/`. Numbers above were collected on Apple Silicon macOS with `amatrix.mlx` installed.

### Flagship 1: Shared-X Many-Y QR (8.6× at 128 RHS)

```sh
Rscript tools/benchmark-flagship-many-y.R
```

Vary `rhs_cols` in the script to reproduce the table above. The 8.6× number corresponds to `X = 1024×128`, `rhs_cols = 128`, MLX backend.

### Flagship 2: GPU Cholesky + Batched Solve (3.4× factor, 2.8× end-to-end)

```sh
Rscript -e 'source("tools/benchmark-cholesky-runtime.R", local = TRUE)'
```

Direct `Rscript tools/benchmark-cholesky-runtime.R` may trip an MLX startup crash on some Apple Silicon setups — use the `-e source(...)` form. The script covers ridge-like (`768×768`, 64 RHS) and kernel-like (`640×640`, 32 RHS) SPD workloads.

### Flagship 3: RSVD on Large Matrices

```sh
Rscript -e 'source("tools/benchmark-svd-factor.R", local = TRUE)'
```

For MLX steady-state SVD timing, see the commands documented in
[docs/gpu-svd-analysis.md](gpu-svd-analysis.md). The calibration grid can be
printed with `tools/print-svd-factor-calibration.R`.

### Summary Table

| Workflow | Matrix | RHS | CPU | MLX | Speedup |
|---|---|---|---|---|---|
| QR many-Y    | 1024×128 | 128 | 0.0215 s | 0.0025 s | **8.6×** |
| Cholesky factor | 768×768 | 64 | 0.065 s | 0.019 s | **3.4×** |
| Cholesky solve  | 768×768 | 64 | 0.080 s | 0.029 s | **2.8×** |
| Cholesky factor | 640×640 | 32 | 0.036 s | 0.011 s | **3.3×** |

All numbers are wall-clock single-run medians on Apple Silicon. MLX backend uses `mode = "fast"` (float32). CPU baseline uses base R `qr()` / `chol()` / `backsolve()` with the shared-X cache enabled.
