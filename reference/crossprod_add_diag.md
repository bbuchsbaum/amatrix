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
#> [1,]  3.7282371 -0.8952379  1.9052199 -2.3326109
#> [2,] -0.8952379  2.4858284 -1.9230556  0.6089594
#> [3,]  1.9052199 -1.9230556  5.1833130 -0.6656873
#> [4,] -2.3326109  0.6089594 -0.6656873  8.4926539
```
