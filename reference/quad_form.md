# Evaluate a quadratic form using a Cholesky factor

Computes `t(v) %*% solve(X) %*% v` (for a vector `v`) or
`t(V) %*% solve(X) %*% V` (for a matrix `V`) efficiently via the
Cholesky factor of `X`, without forming the inverse.

## Usage

``` r
quad_form(factor, v)
```

## Arguments

- factor:

  An
  [`amChol`](https://bbuchsbaum.github.io/amatrix/reference/amChol-class.md)
  object from
  [`chol_factor`](https://bbuchsbaum.github.io/amatrix/reference/chol_factor.md).

- v:

  Numeric vector or matrix. For a vector, the result is a scalar; for a
  matrix with `p` columns, the result is a `p x p` matrix.

## Value

Scalar double (when `v` is a vector) or numeric matrix of dimensions
`ncol(v) x ncol(v)` containing the quadratic form.

## Examples

``` r
m <- crossprod(matrix(rnorm(16), 4, 4)) + diag(4)
A <- adgeMatrix(m)
fac <- chol_factor(A)
v <- rnorm(4)
quad_form(fac, v)
#> [1] 2.143237
```
