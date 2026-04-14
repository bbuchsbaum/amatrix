# Cholesky factorization for adgCMatrix

Compute the Cholesky factorization of a sparse symmetric
positive-definite `adgCMatrix`.

## Usage

``` r
# S4 method for class 'adgCMatrix'
chol(x, ...)
```

## Arguments

- x:

  A symmetric positive-definite `adgCMatrix`.

- ...:

  Further arguments passed to
  [`Matrix::chol`](https://rdrr.io/pkg/Matrix/man/chol-methods.html).

## Value

An `adgCMatrix` or sparse Cholesky factor object.
