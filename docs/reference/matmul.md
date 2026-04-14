# Matrix multiplication

Multiplies two matrices, routing to an accelerated backend when
available. Plain numeric vectors supplied as `y` are promoted to a
column matrix and the result is dropped back to a vector.

## Usage

``` r
matmul(x, y)
```

## Arguments

- x:

  A matrix or `aMatrix` object.

- y:

  A matrix, `aMatrix` object, or numeric vector.

## Value

A matrix (or numeric vector when `y` was a vector) of dimensions
`nrow(x)` by `ncol(y)`.

## Examples

``` r
A <- adgeMatrix(matrix(1:6, 2, 3))
B <- adgeMatrix(matrix(1:6, 3, 2))
matmul(A, B)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 2 x 2 Matrix of class "adgeMatrix"
#>      [,1] [,2]
#> [1,]   22   49
#> [2,]   28   64
```
