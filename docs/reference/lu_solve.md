# Solve a linear system using an LU factor

Solves `A %*% x = B` where `A` is the square matrix stored in `factor`,
delegating to [`base::solve`](https://rdrr.io/r/base/solve.html) (LAPACK
`DGESV`).

## Usage

``` r
lu_solve(factor, B)
```

## Arguments

- factor:

  An
  [`amLU`](https://bbuchsbaum.github.io/amatrix/reference/amLU-class.md)
  object from
  [`lu_factor`](https://bbuchsbaum.github.io/amatrix/reference/lu_factor.md).

- B:

  Numeric vector or matrix; the right-hand side. The number of rows must
  equal the dimension of `factor@A`.

## Value

Numeric vector or matrix `x` satisfying `A %*% x == B`. Returns a vector
when `B` is a vector or single-column matrix.

## Examples

``` r
m <- matrix(c(2, 1, 5, 3), nrow = 2)
fac <- lu_factor(m)
b <- c(1, 2)
lu_solve(fac, b)
#> [1] -7  3
```
