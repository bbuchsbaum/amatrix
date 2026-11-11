# LU-factor a general square matrix

Computes the pivoted LU factorization of a square numeric matrix or
`adgeMatrix` (LAPACK `DGETRF`) and returns it in an
[`amLU`](https://bbuchsbaum.github.io/amatrix/reference/amLU-class.md)
object for repeated use by
[`lu_solve`](https://bbuchsbaum.github.io/amatrix/reference/lu_solve.md).

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
