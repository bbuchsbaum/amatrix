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
#>             [,1]        [,2]        [,3]       [,4]       [,5]
#> [1,]  0.64444500  0.16630568  0.08517043 -0.2312366  0.8627601
#> [2,]  0.16630568  0.15443738 -0.02017602 -0.3607475  0.1652764
#> [3,]  0.08517043 -0.02017602  0.17693276 -0.1520681 -0.1489719
#> [4,] -0.23123658 -0.36074745 -0.15206813  1.3279877  0.4241809
#> [5,]  0.86276011  0.16527644 -0.14897187  0.4241809  2.0028797
```
