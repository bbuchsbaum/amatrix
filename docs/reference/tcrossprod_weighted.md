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
#>            [,1]       [,2]       [,3]        [,4]        [,5]
#> [1,]  0.9432699 -1.1333390 -0.7852048  0.21146735 -0.39017332
#> [2,] -1.1333390  1.4618277  0.8465792 -0.18923928  0.52475904
#> [3,] -0.7852048  0.8465792  1.1118229 -0.49358383  0.23267495
#> [4,]  0.2114674 -0.1892393 -0.4935838  0.64795387 -0.08271198
#> [5,] -0.3901733  0.5247590  0.2326749 -0.08271198  0.20548724
```
