# Store a general square matrix for LU-based solving

Wraps a square numeric matrix or `adgeMatrix` in an
[`amLU`](https://bbuchsbaum.github.io/amatrix/reference/amLU-class.md)
object. The actual LU decomposition is performed by
[`base::solve`](https://rdrr.io/r/base/solve.html) at solve time via
LAPACK's `DGESV`.

## Usage

``` r
lu_factor(A)
```

## Arguments

- A:

  A square numeric matrix or `adgeMatrix`.

## Value

An
[`amLU`](https://bbuchsbaum.github.io/amatrix/reference/amLU-class.md)
object.

## Examples

``` r
m <- matrix(c(2, 1, 5, 3), nrow = 2)
fac <- lu_factor(m)
fac
#> amLU [2x2 | NA | source: NA]
```
