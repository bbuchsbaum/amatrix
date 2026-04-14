# Project new data onto SVD left singular vectors

Computes `t(U) %*% Y`, where `U` is the matrix of left singular vectors
stored in `factor`. Routes through a GPU backend when the factor was
computed in `"fast"` precision and a device copy of `t(U)` is available.

## Usage

``` r
svd_project(factor, Y)
```

## Arguments

- factor:

  An
  [`amSVD`](https://bbuchsbaum.github.io/amatrix/reference/amSVD-class.md)
  object from
  [`svd_factor`](https://bbuchsbaum.github.io/amatrix/reference/svd_factor.md).

- Y:

  Numeric matrix or vector with `nrow(Y)` equal to the number of rows of
  the original source matrix.

## Value

Numeric matrix of dimensions `k x ncol(Y)` containing the projected
coordinates.

## Examples

``` r
m <- matrix(rnorm(30), nrow = 6)
A <- adgeMatrix(m)
fac <- svd_factor(A, k = 2L)
coords <- svd_project(fac, m)
```
