# Reconstruct data from SVD coordinates

Computes `V %*% (Z / d)`, mapping coordinates in the rank-`k` SVD
subspace back to the original column space. Routes through a GPU backend
when the factor was computed in `"fast"` precision and a device copy of
`V` is available.

## Usage

``` r
svd_reconstruct(factor, Z)
```

## Arguments

- factor:

  An
  [`amSVD`](https://bbuchsbaum.github.io/amatrix/reference/amSVD-class.md)
  object from
  [`svd_factor`](https://bbuchsbaum.github.io/amatrix/reference/svd_factor.md).

- Z:

  Numeric matrix or vector with `nrow(Z)` equal to `factor@k`.

## Value

Numeric matrix of dimensions `ncol(X) x ncol(Z)`, where `X` is the
original source matrix.

## Examples

``` r
m <- matrix(rnorm(30), nrow = 6)
A <- adgeMatrix(m)
fac <- svd_factor(A, k = 2L)
coords <- svd_project(fac, m)
approx <- svd_reconstruct(fac, coords)
```
