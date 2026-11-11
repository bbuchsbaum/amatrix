# Leave-one-out cross-validation for linear models

Computes exact leave-one-out (LOO) prediction errors from a single QR
factorization, using the identity \\e\_{(i)} = e_i / (1 - h\_{ii})\\
where \\e_i\\ are the ordinary residuals and \\h\_{ii}\\ the leverages.
The result is identical to refitting the model `n` times with one
observation dropped each time. `X` must have full column rank.

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

  Retained for compatibility; only `"qr"` is supported.

- ...:

  Additional arguments forwarded to
  [`qr`](https://rdrr.io/r/base/qr.html).

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
X <- adgeMatrix(matrix(rnorm(50), nrow = 10))
y <- rnorm(10)
cv <- lm_loo_cv(X, y)
cv$mse
#> [1] 0.9893432
```
