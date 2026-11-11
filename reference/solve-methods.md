# Solve a linear system for adgeMatrix

Compute the solution to \\a x = b\\ or the matrix inverse of `a` when
`b` is missing, dispatching through the amatrix backend.

## Usage

``` r
# S4 method for class 'adgeMatrix,missing'
solve(a, b, ...)

# S4 method for class 'adgeMatrix,ANY'
solve(a, b, ...)
```

## Arguments

- a:

  An `adgeMatrix` coefficient matrix.

- b:

  A matrix or vector right-hand side, or missing for matrix inversion.

- ...:

  Further arguments passed to the backend.

## Value

An `adgeMatrix` (or numeric vector when `b` is a plain vector)
containing the solution.

## Examples

``` r
A <- adgeMatrix(crossprod(matrix(rnorm(9), 3, 3)) + 3 * diag(3))
solve(A)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 3 x 3 Matrix of class "adgeMatrix"
#>            [,1]       [,2]       [,3]
#> [1,] 0.17437539 0.01955487 0.06316321
#> [2,] 0.01955487 0.16066089 0.01948869
#> [3,] 0.06316321 0.01948869 0.11838233
```
