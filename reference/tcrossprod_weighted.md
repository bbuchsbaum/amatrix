# Weighted outer cross-product XWX'

Computes \\X \mathrm{diag}(w) X^T\\, an `n x n` weighted outer
cross-product. A GPU-resident fast path is used when available.

## Usage

``` r
tcrossprod_weighted(X, w)
```

## Arguments

- X:

  Numeric matrix or `adgeMatrix` of shape `[n, p]`.

- w:

  Positive numeric vector of length `n`; observation weights.

## Value

An `adgeMatrix` of shape `[n, n]`.

## See also

[`crossprod_weighted`](https://bbuchsbaum.github.io/amatrix/reference/crossprod_weighted.md),
[`xty_weighted`](https://bbuchsbaum.github.io/amatrix/reference/xty_weighted.md)

## Examples

``` r
X <- matrix(rnorm(20), nrow = 5)
w <- runif(5)
tcrossprod_weighted(X, w)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 5 x 5 Matrix of class "adgeMatrix"
#>             [,1]        [,2]       [,3]       [,4]       [,5]
#> [1,]  0.55589566 -0.07856360 -0.4222609 0.01385347  0.2243879
#> [2,] -0.07856360  0.32753137  0.5437550 0.06956964 -0.1420277
#> [3,] -0.42226094  0.54375505  1.4718267 0.40334292 -0.1996990
#> [4,]  0.01385347  0.06956964  0.4033429 0.36659923  0.1036148
#> [5,]  0.22438785 -0.14202768 -0.1996990 0.10361475  0.1795020
```
