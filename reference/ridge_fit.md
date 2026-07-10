# Fit a single ridge regression model

Solves the penalized least-squares problem \\\min\_\beta \\Y -
X\beta\\^2 + \lambda \\\beta\\^2\\ for a single penalty value `lambda`.

## Usage

``` r
ridge_fit(
  X,
  Y,
  lambda,
  intercept = FALSE,
  penalize_intercept = FALSE,
  include_fitted = TRUE,
  include_residuals = TRUE,
  cache = TRUE
)
```

## Arguments

- X:

  Numeric matrix or `adgeMatrix` of predictors, shape `[n, p]`.

- Y:

  Numeric matrix, vector, or `adgeMatrix` of responses, shape `[n, q]`.

- lambda:

  Non-negative scalar ridge penalty.

- intercept:

  Logical; when `TRUE` a column of ones is prepended to `X` before
  fitting.

- penalize_intercept:

  Logical; when `FALSE` (default) the intercept coefficient is excluded
  from the penalty.

- include_fitted:

  Logical; include fitted values in the result.

- include_residuals:

  Logical; include residuals in the result.

- cache:

  Logical; cache \\X^T X\\ for reuse across calls sharing the same `X`.

## Value

An object of class `"ridge_fit"`, a named list containing:

- coefficients:

  `adgeMatrix` of shape `[p, q]`.

- fitted.values:

  `adgeMatrix` of shape `[n, q]`, or `NULL` when
  `include_fitted = FALSE`.

- residuals:

  `adgeMatrix` of shape `[n, q]`, or `NULL` when
  `include_residuals = FALSE`.

- lambda:

  The penalty value used.

- rank:

  Integer model rank.

- df.residual:

  Residual degrees of freedom.

## See also

[`ridge_path`](https://bbuchsbaum.github.io/amatrix/reference/ridge_path.md),
[`lm_fit`](https://bbuchsbaum.github.io/amatrix/reference/lm_fit.md)

## Examples

``` r
X <- matrix(rnorm(50), nrow = 10)
y <- rnorm(10)
fit <- ridge_fit(X, y, lambda = 1)
coef(fit)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 5 x 1 Matrix of class "adgeMatrix"
#>             [,1]
#> [1,] -0.18779241
#> [2,] -0.36505283
#> [3,] -0.08561205
#> [4,] -0.17422701
#> [5,]  0.17126016
```
