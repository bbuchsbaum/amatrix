# Leave-one-out cross-validation for linear models

Computes exact leave-one-out (LOO) prediction errors by refitting the
model `n` times, each time dropping one observation. Uses `qr_downdate`
internally to avoid recomputing the full factorization from scratch.

## Usage

``` r
lm_loo_cv(X, y, method = "qr", ...)
```

## Arguments

- X:

  Numeric matrix or `adgeMatrix` of predictors, shape `[n, p]`. No
  intercept column is added automatically.

- y:

  Numeric vector or single-column matrix of responses, length `n`.

- method:

  Character string passed to `am_qr` controlling the QR algorithm.
  Currently only `"qr"` is supported.

- ...:

  Additional arguments forwarded to `am_qr`.

## Value

A named list with two elements:

- residuals:

  Numeric vector of length `n`: LOO prediction errors \\y_i -
  \hat{y}\_i^{(-i)}\\.

- mse:

  Scalar mean squared LOO error.

## See also

[`lm_fit`](https://bbuchsbaum.github.io/amatrix/reference/lm_fit.md),
[`many_lm`](https://bbuchsbaum.github.io/amatrix/reference/many_lm.md)

## Examples

``` r
if (FALSE) { # \dontrun{
X <- adgeMatrix(matrix(rnorm(50), nrow = 10))
y <- rnorm(10)
cv <- lm_loo_cv(X, y)
cv$mse
} # }
```
