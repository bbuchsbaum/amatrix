# Weighted cross-product X'Wy

Computes \\X^T \mathrm{diag}(w) y\\, a `p x k` weighted cross-product
between `X` and response matrix `y`. A GPU-resident fast path is used
when available.

## Usage

``` r
xty_weighted(X, w, y)
```

## Arguments

- X:

  Numeric matrix or `adgeMatrix` of shape `[n, p]`.

- w:

  Positive numeric vector of length `n`; observation weights.

- y:

  Numeric vector or matrix of shape `[n, k]`; response(s).

## Value

An `adgeMatrix` of shape `[p, k]`.

## See also

[`crossprod_weighted`](https://bbuchsbaum.github.io/amatrix/reference/crossprod_weighted.md),
[`tcrossprod_weighted`](https://bbuchsbaum.github.io/amatrix/reference/tcrossprod_weighted.md)

## Examples

``` r
X <- matrix(rnorm(20), nrow = 5)
w <- runif(5)
y <- rnorm(5)
xty_weighted(X, w, y)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 4 x 1 Matrix of class "adgeMatrix"
#>              [,1]
#> [1,]  1.209024440
#> [2,] -0.001972357
#> [3,] -1.227735079
#> [4,]  1.016346260
```
