# Compute a ridge regression solution path

Fits ridge regression for every penalty value in `lambdas` via a single
thin SVD of `X`, returning coefficients for all penalties at once.

## Usage

``` r
ridge_path(X, Y, lambdas, k = NULL, ...)
```

## Arguments

- X:

  Numeric matrix or `adgeMatrix` of predictors, shape `[n, p]`.

- Y:

  Numeric matrix, vector, or `adgeMatrix` of responses, shape `[n, q]`.

- lambdas:

  Positive numeric vector of ridge penalty values. Must satisfy
  `all(lambdas > 0)`.

- k:

  Integer or `NULL`. Number of singular values to retain in the
  truncated SVD. When `NULL`, defaults to `min(nrow(X), ncol(X))`.

- ...:

  Additional arguments forwarded to `svd_factor`.

## Value

An object of class `"ridge_path"`, a named list containing:

- coef:

  Numeric array of shape `[p, q, length(lambdas)]`; coefficient matrix
  for each penalty.

- lambdas:

  The input penalty vector.

- svd:

  The `amSVD` factor object used internally.

- k:

  Integer number of singular values actually used.

## See also

[`ridge_fit`](https://bbuchsbaum.github.io/amatrix/reference/ridge_fit.md),
[`svd_factor`](https://bbuchsbaum.github.io/amatrix/reference/svd_factor.md)

## Examples

``` r
X <- matrix(rnorm(60), nrow = 15)
y <- rnorm(15)
path <- ridge_path(X, y, lambdas = c(0.1, 1, 10))
dim(path$coef)
#> [1] 4 1 3
```
