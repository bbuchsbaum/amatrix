# Eager Kronecker product

Computes the Kronecker product \\A \otimes B\\ and returns the result as
an `adgeMatrix`. Accepts plain matrices or any `aMatrix` subclass. For a
lazy variant that avoids forming the full product see
[`kron_matrix`](https://bbuchsbaum.github.io/amatrix/reference/kron_matrix.md).

## Usage

``` r
kron(A, B)
```

## Arguments

- A, B:

  Matrices or `aMatrix` objects.

## Value

An `adgeMatrix` of dimension `(nrow(A)*nrow(B)) x (ncol(A)*ncol(B))`.

## See also

[`kron_matrix`](https://bbuchsbaum.github.io/amatrix/reference/kron_matrix.md)

## Examples

``` r
A <- matrix(1:4, 2, 2)
B <- diag(2)
kron(A, B)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 4 x 4 Matrix of class "adgeMatrix"
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    0    3    0
#> [2,]    0    1    0    3
#> [3,]    2    0    4    0
#> [4,]    0    2    0    4
```
