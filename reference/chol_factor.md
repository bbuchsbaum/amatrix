# Compute the Cholesky factorization of an adgeMatrix

Computes the upper-triangular Cholesky factor `R` of a symmetric
positive-definite `adgeMatrix` `X` such that `t(R) %*% R == X`. Results
are cached by `object_id`; repeated calls with the same object return
the cached factor.

## Usage

``` r
chol_factor(X)
```

## Arguments

- X:

  An `adgeMatrix` that is symmetric positive definite.

## Value

An
[`amChol`](https://bbuchsbaum.github.io/amatrix/reference/amChol-class.md)
object.

## Examples

``` r
m <- crossprod(matrix(rnorm(16), 4, 4)) + diag(4)
A <- adgeMatrix(m)
fac <- chol_factor(A)
fac
#> amChol [4x4 | strict | source: 20260710061154.331847-58383fa0:am:30]
```
