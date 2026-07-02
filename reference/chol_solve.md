# Solve a linear system using a Cholesky factor

Solves `X %*% x = B` where `X` is the symmetric positive-definite matrix
whose Cholesky factorization is stored in `factor`. Dispatches to a GPU
backend when the factor was computed in `"fast"` precision and a
device-resident factor is available.

## Usage

``` r
chol_solve(factor, B)
```

## Arguments

- factor:

  An
  [`amChol`](https://bbuchsbaum.github.io/amatrix/reference/amChol-class.md)
  object from
  [`chol_factor`](https://bbuchsbaum.github.io/amatrix/reference/chol_factor.md).

- B:

  Numeric vector or matrix; the right-hand side. The number of rows must
  equal the dimension of the factor.

## Value

Numeric vector or matrix `x` satisfying `X %*% x == B`. Returns a vector
when `B` is a vector.

## Examples

``` r
m <- crossprod(matrix(rnorm(16), 4, 4)) + diag(4)
A <- adgeMatrix(m)
fac <- chol_factor(A)
b <- rnorm(4)
x <- chol_solve(fac, b)
```
