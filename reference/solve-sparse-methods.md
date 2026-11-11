# Solve a linear system for adgCMatrix

Compute the solution to \\a x = b\\ or the inverse of `a` when `b` is
missing, for a sparse `adgCMatrix` coefficient matrix.

## Usage

``` r
# S4 method for class 'adgCMatrix,missing'
solve(a, b, ...)

# S4 method for class 'adgCMatrix,ANY'
solve(a, b, ...)
```

## Arguments

- a:

  An `adgCMatrix` coefficient matrix.

- b:

  A matrix or vector right-hand side, or missing for inversion.

- ...:

  Further arguments passed to the backend.

## Value

An `adgeMatrix` or `adgCMatrix` containing the solution.

## Examples

``` r
S <- Matrix::sparseMatrix(i = c(1, 2, 3, 1), j = c(1, 2, 3, 3),
                          x = c(2, 3, 4, 1))
A <- adgCMatrix(S)
solve(A, c(1, 2, 3))
#> [1] 0.1250000 0.6666667 0.7500000
```
