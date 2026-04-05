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

- `rhs_cols = 8`
  - base cached QR: about `0.0016 s`
  - MLX native resident: about `0.0020 s`
- `rhs_cols = 32`
  - base cached QR: about `0.0058 s`
  - MLX native resident: about `0.0018 s`
- `rhs_cols = 128`
  - base cached QR: about `0.0215 s`
  - MLX native resident: about `0.0025 s`

That is the clearest current reason to adopt the model surface.

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
| `"balanced"` | strict float64, auto routing — GPU where numerically safe (default) |
| `"fast"` | float32-oriented, auto routing — full GPU throughput |

Default (`mode` omitted) uses `"balanced"` semantics: strict precision, CPU unless a backend is specified.

For accelerator-oriented execution, specify both `mode` and `backend`:

```r
X <- adgeMatrix(design, mode = "fast", backend = "mlx")
fit <- am_many_lm(X, Y_many, method = "qr", cache = TRUE)
```

The `backend=` argument is an escape hatch for users who know which accelerator they want. Without it, the system stays on CPU regardless of `mode`.

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
