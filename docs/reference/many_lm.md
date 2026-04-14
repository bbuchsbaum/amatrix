# Fit multiple linear models against a shared design matrix

The flagship batch regression function. Solves \\Y \approx X B\\ where
`Y` is a matrix whose columns are independent response variables. When
no GPU backend is active the normal-equations or QR path runs on CPU.
With an active GPU backend the QR path dispatches the factorization to
the device and keeps intermediate results resident to minimise host
round-trips.

## Usage

``` r
many_lm(
  X,
  Y,
  weights = NULL,
  intercept = FALSE,
  include_fitted = FALSE,
  include_residuals = FALSE,
  cache = TRUE,
  method = c("normal", "qr")
)
```

## Arguments

- X:

  Numeric matrix or `adgeMatrix` of predictors, shape `[n, p]`.

- Y:

  Numeric matrix or `adgeMatrix` of responses, shape `[n, q]`. Each
  column is fitted independently against `X`.

- weights:

  Optional numeric vector of length `n` with non-negative observation
  weights. When supplied, weighted least squares is used.

- intercept:

  Logical. When `TRUE`, a column of ones is prepended to `X` before
  fitting.

- include_fitted:

  Logical. When `TRUE`, fitted values are stored in the returned object.

- include_residuals:

  Logical. When `TRUE`, residuals are stored in the returned object.

- cache:

  Logical. When `TRUE`, the design-matrix factorization is cached for
  reuse when `X` is the same across calls.

- method:

  Solver: `"normal"` (normal equations) or `"qr"` (QR decomposition).
  Ignored when `weights` is non-`NULL` and the weighted path is
  selected.

## Value

An object of class `"am_many_lm_fit"`, a named list containing:

- coefficients:

  `adgeMatrix` of shape `[p, q]`.

- fitted.values:

  `adgeMatrix` `[n, q]` or `NULL`.

- residuals:

  `adgeMatrix` `[n, q]` or `NULL`.

- rss:

  Numeric vector of length `q`: residual sums of squares.

- sigma2:

  Numeric vector of length `q`: residual variances.

- rank:

  Integer model rank.

- df.residual:

  Residual degrees of freedom.

## See also

[`lm_fit`](https://bbuchsbaum.github.io/amatrix/reference/lm_fit.md),
[`array_lm`](https://bbuchsbaum.github.io/amatrix/reference/array_lm.md)

## Examples

``` r
X <- matrix(rnorm(100), nrow = 20)
Y <- matrix(rnorm(60), nrow = 20)
fit <- many_lm(X, Y, method = "qr")
dim(coef(fit))
#> [1] 5 3
```
