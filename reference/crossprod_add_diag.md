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
#>            [,1]       [,2]       [,3]       [,4]
#> [1,]  2.0350580 -0.7507005 -0.2557977  0.2590892
#> [2,] -0.7507005  5.8658975  1.0719273 -2.9140488
#> [3,] -0.2557977  1.0719273  0.8326541 -0.2257762
#> [4,]  0.2590892 -2.9140488 -0.2257762  4.2722464
```
