# Weighted cross-product X'WX

Computes \\X^T \mathrm{diag}(w) X\\, a `p x p` weighted cross-product. A
GPU-resident fast path is used when available.

## Usage

``` r
crossprod_weighted(X, w)
```

## Arguments

- X:

  Numeric matrix or `adgeMatrix` of shape `[n, p]`.

- w:

  Positive numeric vector of length `n`; observation weights.

## Value

An `adgeMatrix` of shape `[p, p]`.

## See also

[`tcrossprod_weighted`](https://bbuchsbaum.github.io/amatrix/reference/tcrossprod_weighted.md),
[`xty_weighted`](https://bbuchsbaum.github.io/amatrix/reference/xty_weighted.md)

## Examples

``` r
X <- matrix(rnorm(20), nrow = 5)
w <- runif(5)
crossprod_weighted(X, w)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 4 x 4 Matrix of class "adgeMatrix"
#>            [,1]       [,2]       [,3]       [,4]
#> [1,]  4.9969074 -2.0771600  0.7519243 -1.6156023
#> [2,] -2.0771600  5.3798637 -0.1886434  3.5719034
#> [3,]  0.7519243 -0.1886434  0.4495651  0.4422088
#> [4,] -1.6156023  3.5719034  0.4422088  6.6543849
```
