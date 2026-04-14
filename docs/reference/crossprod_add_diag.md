# Cross-product plus diagonal perturbation

Computes \\X^T X + \lambda I\\ (scalar `lambda`) or \\X^T X +
\mathrm{diag}(\lambda)\\ (vector `lambda`) in a single fused call.

## Usage

``` r
crossprod_add_diag(X, lambda)
```

## Arguments

- X:

  Numeric matrix or `adgeMatrix` of shape `[n, p]`.

- lambda:

  Scalar or numeric vector of length `p`; diagonal perturbation to add
  to the cross-product.

## Value

An `adgeMatrix` of shape `[p, p]`: the perturbed cross-product.

## See also

[`crossprod_weighted`](https://bbuchsbaum.github.io/amatrix/reference/crossprod_weighted.md)

## Examples

``` r
X <- matrix(rnorm(20), nrow = 5)
crossprod_add_diag(X, lambda = 0.1)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 4 x 4 Matrix of class "adgeMatrix"
#>            [,1]       [,2]      [,3]      [,4]
#> [1,]  2.0673925  0.8072824 -1.911997 -1.522897
#> [2,]  0.8072824  2.7185409 -2.737164  2.760107
#> [3,] -1.9119967 -2.7371637  6.664585 -4.864255
#> [4,] -1.5228973  2.7601068 -4.864255 12.335061
```
