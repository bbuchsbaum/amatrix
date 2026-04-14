# Extract the diagonal of a Cholesky factor

Returns the diagonal of the upper-triangular matrix `R` stored in an
[`amChol`](https://bbuchsbaum.github.io/amatrix/reference/amChol-class.md)
object.

## Usage

``` r
chol_diag(factor)
```

## Arguments

- factor:

  An
  [`amChol`](https://bbuchsbaum.github.io/amatrix/reference/amChol-class.md)
  object from
  [`chol_factor`](https://bbuchsbaum.github.io/amatrix/reference/chol_factor.md).

## Value

Numeric vector of length equal to the matrix dimension.

## Examples

``` r
m <- crossprod(matrix(rnorm(16), 4, 4)) + diag(4)
A <- adgeMatrix(m)
fac <- chol_factor(A)
chol_diag(fac)
#> [1] 2.998567 1.344080 1.722819 1.905674
```
