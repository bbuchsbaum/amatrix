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
