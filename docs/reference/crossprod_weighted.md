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
#> [1,]  1.8667988  1.3585042  0.9274433 -0.7039125
#> [2,]  1.3585042  1.2321113  0.7798011 -0.6539223
#> [3,]  0.9274433  0.7798011  0.5117053 -0.4308621
#> [4,] -0.7039125 -0.6539223 -0.4308621  0.4368487
```
