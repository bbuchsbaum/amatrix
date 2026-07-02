# Row and column diagonal scaling

Scale each row or column of a matrix by a numeric vector, equivalent to
left- or right-multiplying by a diagonal matrix. `rowscale` computes
`diag(d) %*% X` (row \\i\\ scaled by `d[i]`); `colscale` computes
`X %*% diag(d)` (column \\j\\ scaled by `d[j]`).

## Usage

``` r
rowscale(X, d)

colscale(X, d)
```

## Arguments

- X:

  A matrix or `aMatrix` object.

- d:

  Numeric vector of scale factors. Length must equal `nrow(X)` for
  `rowscale` and `ncol(X)` for `colscale`.

## Value

A matrix or `aMatrix` of the same dimensions as `X`.

## Examples

``` r
m <- matrix(1:6, 2, 3)
rowscale(m, c(2, 0.5))
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 2 x 3 Matrix of class "adgeMatrix"
#>      [,1] [,2] [,3]
#> [1,]    2    6   10
#> [2,]    1    2    3
colscale(m, c(1, 2, 3))
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 2 x 3 Matrix of class "adgeMatrix"
#>      [,1] [,2] [,3]
#> [1,]    1    6   15
#> [2,]    2    8   18
```
