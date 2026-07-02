# Cholesky factorization for adgeMatrix

Compute the Cholesky factor of a symmetric positive-definite
`adgeMatrix`, dispatching through the amatrix backend.

## Usage

``` r
# S4 method for class 'adgeMatrix'
chol(x, ...)
```

## Arguments

- x:

  A symmetric positive-definite `adgeMatrix`.

- ...:

  Further arguments passed to the backend.

## Value

An `adgeMatrix` containing the upper triangular Cholesky factor.

## Examples

``` r
S <- adgeMatrix(crossprod(matrix(rnorm(9), 3, 3)) + 3 * diag(3))
R <- chol(S)
```
