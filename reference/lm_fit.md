# Fit a single linear model

Solves the ordinary least-squares problem \\Y \approx X \beta\\ for a
single design matrix `X` and one or more response columns `Y`.

## Usage

``` r
lm_fit(
  X,
  Y,
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

- intercept:

  Logical. When `TRUE`, a column of ones is prepended to `X` before
  fitting.

- include_fitted:

  Logical. When `TRUE`, fitted values are included in the returned
  object.

- include_residuals:

  Logical. When `TRUE`, residuals are included in the returned object.

- cache:

  Logical. When `TRUE`, the \\X^T X\\ or QR factorization is cached for
  reuse across calls sharing the same `X`.

- method:

  Solver method: `"normal"` (normal equations, default) or `"qr"` (QR
  decomposition). The normal equations are the fastest route and reuse a
  cached \\X^\top X\\, but they square the condition number of `X`: a
  warning of class `amatrix_ill_conditioned` is raised when \\X^\top X\\
  is close to singular, in which case `"qr"` should be used. A
  rank-deficient `X` is an error of class `amatrix_rank_deficient`.

## Value

An object of class `"lm_fit"`, a named list containing:

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

[`many_lm`](https://bbuchsbaum.github.io/amatrix/reference/many_lm.md),
[`array_lm`](https://bbuchsbaum.github.io/amatrix/reference/array_lm.md)

## Examples

``` r
X <- matrix(rnorm(50), nrow = 10)
y <- rnorm(10)
fit <- lm_fit(X, y)
coef(fit)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 5 x 1 Matrix of class "adgeMatrix"
#>             [,1]
#> [1,] -0.01768660
#> [2,]  0.54747493
#> [3,] -0.02178289
#> [4,]  0.16732455
#> [5,] -0.09059265
```
