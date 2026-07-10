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
#>             [,1]        [,2]        [,3]        [,4]        [,5]
#> [1,]  0.08939343  0.05945836 -0.07658731 -0.08858679 -0.23780447
#> [2,]  0.05945836  0.89731015 -0.42414500  0.04972706 -0.68340116
#> [3,] -0.07658731 -0.42414500  0.93535717  1.35793648 -0.03207137
#> [4,] -0.08858679  0.04972706  1.35793648  2.66995100 -0.63481289
#> [5,] -0.23780447 -0.68340116 -0.03207137 -0.63481289  1.32606690
```
