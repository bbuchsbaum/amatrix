# amatrix

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/bbuchsbaum/amatrix/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bbuchsbaum/amatrix/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`amatrix` keeps matrix-heavy R code in the `Matrix` idiom while adding backend-aware dispatch and optional accelerator execution. The main target is repeated linear algebra on the same design matrix, especially many-response regression where factorization reuse matters.

## What problem does it solve?

If your workflow already starts with matrices and ends in `%*%`, `crossprod()`, `qr()`, or repeated least-squares fits, `amatrix` gives you:

- Matrix-compatible dense and sparse objects
- predictable CPU fallback when no accelerator backend is available
- optional fast backends such as MLX on Apple Silicon
- shared-design workflows such as `many_lm()`

You do not need to adopt a new tensor framework or rewrite the analysis around a separate object model.

## Quick start

```r
library(amatrix)

set.seed(1)
X <- matrix(rnorm(120 * 6), nrow = 120, ncol = 6)
Y <- matrix(rnorm(120 * 8), nrow = 120, ncol = 8)

# Safe default: strict CPU semantics
X_am <- adgeMatrix(X)
fit <- many_lm(X_am, Y, method = "qr", cache = TRUE)

dim(coef(fit))
fit$sigma2[1:3]
```

On a machine with an available accelerator backend, the workflow stays the same and only the constructor metadata changes.

```r
X_fast <- adgeMatrix(X, mode = "fast", backend = "mlx")
fit_fast <- many_lm(X_fast, Y, method = "qr", cache = TRUE)
```

## How does backend selection work?

`amatrix` separates the matrix object from the execution backend. You can inspect what is available and what will be chosen for a specific operation.

```r
amatrix_backend_status()
amatrix_backend_plan(X_am, "qr")
amatrix_explain(X_am, "qr")
```

If no accelerator backend is available, the same code still runs on CPU.

## When is it worth using?

- one design matrix and many response columns
- repeated fits where QR or other factors can be reused
- workflows that need Matrix-compatible objects instead of a separate tensor API
- Apple Silicon systems where MLX can accelerate the hot path

## Modes

- `adgeMatrix(x)` uses the conservative default path and is safe on CPU-only machines.
- `adgeMatrix(x, mode = "exact")` pins execution to strict CPU semantics.
- `adgeMatrix(x, mode = "fast", backend = "mlx")` allows reduced-precision accelerator execution when the backend supports it.

## Installation

```r
# Core package
pak::pkg_install("amatrix")

# Optional MLX backend for Apple Silicon
pak::pkg_install(c("amatrix", "amatrix.mlx"))
```

Other accelerator backends can be installed separately when needed.

## Start here

- `vignette("amatrix")` for the getting-started workflow
- `?adgeMatrix` for matrix constructors
- `?many_lm` for the flagship shared-design regression path
- `?amatrix_backend_status` for backend availability and capability checks
