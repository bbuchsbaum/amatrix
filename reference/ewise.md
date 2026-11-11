# Element-wise operations

Apply an element-wise arithmetic operation to one or two matrices,
dispatching to the preferred GPU backend when available.

## Usage

``` r
ewise(op, e1, e2 = NULL)
```

## Arguments

- op:

  Character string: `"+"`, `"-"`, `"*"`, or `"/"`.

- e1:

  A numeric matrix or `adgeMatrix`.

- e2:

  A numeric matrix, `adgeMatrix`, or `NULL` for unary ops.

## Value

A matrix of the same dimensions as `e1`.

## Examples

``` r
x <- adgeMatrix(matrix(1:6, 2, 3))
ewise("+", x, 1)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 2 x 3 Matrix of class "adgeMatrix"
#>      [,1] [,2] [,3]
#> [1,]    2    4    6
#> [2,]    3    5    7
ewise("*", x, x)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 2 x 3 Matrix of class "adgeMatrix"
#>      [,1] [,2] [,3]
#> [1,]    1    9   25
#> [2,]    4   16   36
ewise("exp", x)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 2 x 3 Matrix of class "adgeMatrix"
#>          [,1]     [,2]     [,3]
#> [1,] 2.718282 20.08554 148.4132
#> [2,] 7.389056 54.59815 403.4288
```
