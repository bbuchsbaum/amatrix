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
#>           [,1]      [,2]       [,3]       [,4]
#> [1,] 2.9453675 1.5498585  1.1536397  0.5179425
#> [2,] 1.5498585 1.2960658  0.1666179  0.7701780
#> [3,] 1.1536397 0.1666179  1.7368843 -0.5174326
#> [4,] 0.5179425 0.7701780 -0.5174326  1.3451171
```
