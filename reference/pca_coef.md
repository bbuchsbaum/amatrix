# Principal-component regression coefficients from a truncated SVD

Convenience wrapper that calls
[`svd_project`](https://bbuchsbaum.github.io/amatrix/reference/svd_project.md)
followed by
[`svd_reconstruct`](https://bbuchsbaum.github.io/amatrix/reference/svd_reconstruct.md),
i.e. computes \\V D^{-1} U^\top Y\\: the minimum-norm least-squares
coefficients of `Y` on the rank-`k` approximation of the original matrix
`X`. Multiply `X` by the result to obtain fitted values. Singular values
that are zero to working precision are dropped (pseudo-inverse).

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

Numeric matrix of dimensions `ncol(X) x ncol(Y)`: the coefficient
matrix.

## Examples

``` r
m <- matrix(rnorm(30), nrow = 6)
A <- adgeMatrix(m)
fac <- svd_factor(A, k = 2L)
y <- rnorm(6)
beta <- pca_coef(fac, y)
fitted <- m %*% beta
```
