# Fit linear models with array-shaped response

Wraps
[`many_lm`](https://bbuchsbaum.github.io/amatrix/reference/many_lm.md)
to accept a response `Y` with more than two dimensions (e.g. a 3-D
array). The trailing dimensions are collapsed into columns for fitting
and optionally restored in the output.

## Usage

``` r
array_lm(
  X,
  Y,
  weights = NULL,
  intercept = FALSE,
  include_fitted = FALSE,
  include_residuals = FALSE,
  cache = TRUE,
  method = c("normal", "qr"),
  restore_array = TRUE
)
```

## Arguments

- X:

  Numeric matrix or `adgeMatrix` of predictors, shape `[n, p]`.

- Y:

  Numeric array or matrix of responses. The first dimension must equal
  `n` (observations). Additional dimensions are treated as independent
  response variables.

- weights:

  Optional numeric vector of length `n` with non-negative observation
  weights.

- intercept:

  Logical. When `TRUE`, a column of ones is prepended to `X` before
  fitting.

- include_fitted:

  Logical. When `TRUE`, fitted values are stored in the returned object.

- include_residuals:

  Logical. When `TRUE`, residuals are stored in the returned object.

- cache:

  Logical. When `TRUE`, the design-matrix factorization is cached for
  reuse.

- method:

  Solver: `"normal"` or `"qr"`.

- restore_array:

  Logical. When `TRUE` (default), `rss`, `sigma2`, fitted values, and
  residuals are reshaped to match the original trailing dimensions of
  `Y`.

## Value

An object of class `"am_array_lm_fit"`, a named list containing the same
fields as
[`many_lm`](https://bbuchsbaum.github.io/amatrix/reference/many_lm.md)
plus:

- response_dims:

  Integer vector of trailing dimensions of `Y`.

- rss:

  Array or vector of residual sums of squares.

- sigma2:

  Array or vector of residual variances.

## See also

[`many_lm`](https://bbuchsbaum.github.io/amatrix/reference/many_lm.md),
[`lm_fit`](https://bbuchsbaum.github.io/amatrix/reference/lm_fit.md)

## Examples

``` r
X <- matrix(rnorm(50), nrow = 10)
Y <- array(rnorm(10 * 3 * 4), dim = c(10, 3, 4))
fit <- array_lm(X, Y)
dim(fit$sigma2)
#> NULL
```
