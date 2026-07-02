# Fit a weighted least squares model

Solves the weighted least-squares problem \\\min\_\beta \sum_i w_i
(y_i - x_i^T \beta)^2\\ using either the normal equations or a QR
decomposition of the weight-scaled design.

## Usage

``` r
wls_fit(
  X,
  Y,
  weights,
  intercept = FALSE,
  include_fitted = TRUE,
  include_residuals = TRUE,
  cache = TRUE,
  method = c("normal", "qr")
)
```

## Arguments

- X:

  Numeric matrix or `adgeMatrix` of predictors, shape `[n, p]`.

- Y:

  Numeric matrix, vector, or `adgeMatrix` of responses, shape `[n, q]`.

- weights:

  Positive numeric vector of length `n`; observation weights.

- intercept:

  Logical; when `TRUE` a column of ones is prepended to `X` before
  fitting.

- include_fitted:

  Logical; include fitted values in the result.

- include_residuals:

  Logical; include residuals in the result.

- cache:

  Logical; cache intermediate factorizations for reuse across calls
  sharing the same `X` and `weights`.

- method:

  Solver method: `"normal"` (weighted normal equations, default) or
  `"qr"` (QR on the weight-scaled design).

## Value

An object of class `"wls_fit"`, a named list containing:

- coefficients:

  `adgeMatrix` of shape `[p, q]`.

- fitted.values:

  `adgeMatrix` of shape `[n, q]`, or `NULL` when
  `include_fitted = FALSE`.

- residuals:

  `adgeMatrix` of shape `[n, q]`, or `NULL` when
  `include_residuals = FALSE`.

- rank:

  Integer model rank.

- df.residual:

  Residual degrees of freedom.

## See also

[`lm_fit`](https://bbuchsbaum.github.io/amatrix/reference/lm_fit.md),
[`crossprod_weighted`](https://bbuchsbaum.github.io/amatrix/reference/crossprod_weighted.md)

## Examples

``` r
X <- matrix(rnorm(50), nrow = 10)
y <- rnorm(10)
w <- runif(10, 0.5, 2)
fit <- wls_fit(X, y, weights = w)
coef(fit)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 5 x 1 Matrix of class "adgeMatrix"
#>             [,1]
#> [1,]  0.27670490
#> [2,] -0.32151899
#> [3,] -0.04981365
#> [4,] -0.13500564
#> [5,] -0.74870841
```
