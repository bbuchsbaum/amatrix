# Construct a lazy Kronecker product

Creates a
[`KronMatrix`](https://bbuchsbaum.github.io/amatrix/reference/KronMatrix-class.md)
that stores the factor matrices `A` and `B` without materializing the
full Kronecker product. Standard operations such as `%*%`, `crossprod`,
`solve`, and `determinant` are available and exploit the Kronecker
structure.

## Usage

``` r
kron_matrix(A, B)
```

## Arguments

- A:

  Numeric matrix or object coercible via
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html); the left
  Kronecker factor.

- B:

  Numeric matrix or object coercible via
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html); the right
  Kronecker factor.

## Value

A
[`KronMatrix`](https://bbuchsbaum.github.io/amatrix/reference/KronMatrix-class.md)
of implicit dimensions `c(nrow(A) * nrow(B), ncol(A) * ncol(B))`.

## Examples

``` r
A <- matrix(1:4, 2, 2)
B <- diag(3)
K <- kron_matrix(A, B)
dim(K)
#> [1] 6 6
as.matrix(K)
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    1    0    0    3    0    0
#> [2,]    0    1    0    0    3    0
#> [3,]    0    0    1    0    0    3
#> [4,]    2    0    0    4    0    0
#> [5,]    0    2    0    0    4    0
#> [6,]    0    0    2    0    0    4
```
