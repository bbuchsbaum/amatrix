# Covariance-to-correlation methods for amatrix objects

Bridge [`cov2cor()`](https://rdrr.io/r/stats/cor.html) through Matrix's
covariance-to-correlation methods so standard workflows such as
`cov2cor(crossprod(X))` keep working when
[`crossprod()`](https://rdrr.io/r/base/crossprod.html) preserves an
amatrix class.

## Usage

``` r
# S4 method for class 'adgeMatrix'
cov2cor(V)

# S4 method for class 'adgCMatrix'
cov2cor(V)
```

## Arguments

- V:

  A square `adgeMatrix` or `adgCMatrix`.

## Value

A base R correlation matrix, matching
[`stats::cov2cor()`](https://rdrr.io/r/stats/cor.html) on the
corresponding host matrix.

## Examples

``` r
X <- adgeMatrix(matrix(1:9 + 0, 3, 3))
cov2cor(crossprod(X))
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.9746318 0.9594119
#> [2,] 0.9746318 1.0000000 0.9981909
#> [3,] 0.9594119 0.9981909 1.0000000
```
