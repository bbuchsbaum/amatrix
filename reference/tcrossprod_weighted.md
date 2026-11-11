# Row-weighted outer cross-product

Computes \\W^{1/2} X X^T W^{1/2}\\ with \\W = \mathrm{diag}(w)\\: the
`n x n` outer cross-product of the row-weighted matrix \\W^{1/2} X\\. A
GPU-resident fast path is used when available.

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
#>            [,1]         [,2]       [,3]       [,4]         [,5]
#> [1,]  0.5528871 -0.465195025 -0.0662382 -0.9561148 -0.271782917
#> [2,] -0.4651950  0.903362448  0.1149791  0.5045299 -0.008577877
#> [3,] -0.0662382  0.114979103  0.9767493 -0.7482508 -1.421692970
#> [4,] -0.9561148  0.504529927 -0.7482508  2.6936906  1.982184869
#> [5,] -0.2717829 -0.008577877 -1.4216930  1.9821849  2.498311341
```
