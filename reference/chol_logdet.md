# Log-determinant from a Cholesky factor

Computes `log(det(X))` from the Cholesky factor of a symmetric
positive-definite matrix `X` as `2 * sum(log(diag(R)))`, which avoids
forming the full determinant and is numerically stable.

## Usage

``` r
chol_logdet(factor)
```

## Arguments

- factor:

  An
  [`amChol`](https://bbuchsbaum.github.io/amatrix/reference/amChol-class.md)
  object from
  [`chol_factor`](https://bbuchsbaum.github.io/amatrix/reference/chol_factor.md).

## Value

Scalar double; the log-determinant of the source matrix.

## Examples

``` r
m <- crossprod(matrix(rnorm(16), 4, 4)) + diag(4)
A <- adgeMatrix(m)
fac <- chol_factor(A)
chol_logdet(fac)
#> [1] 5.853034
```
