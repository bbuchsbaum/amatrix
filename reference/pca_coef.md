# Project and reconstruct data using a truncated SVD

Convenience wrapper that calls
[`svd_project`](https://bbuchsbaum.github.io/amatrix/reference/svd_project.md)
followed by
[`svd_reconstruct`](https://bbuchsbaum.github.io/amatrix/reference/svd_reconstruct.md),
yielding the rank-`k` least-squares approximation of `Y` in the column
space of the original matrix.

## Usage

``` r
pca_coef(factor, Y)
```

## Arguments

- factor:

  An
  [`amSVD`](https://bbuchsbaum.github.io/amatrix/reference/amSVD-class.md)
  object from
  [`svd_factor`](https://bbuchsbaum.github.io/amatrix/reference/svd_factor.md).

- Y:

  Numeric matrix or vector to project and reconstruct. Must have
  `nrow(Y)` equal to the number of rows of the original source matrix.

## Value

Numeric matrix with the same dimensions as `Y`, containing the rank-`k`
approximation.

## Examples

``` r
m <- matrix(rnorm(30), nrow = 6)
A <- adgeMatrix(m)
fac <- svd_factor(A, k = 2L)
approx <- pca_coef(fac, m)
```
