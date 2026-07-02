# Matrix multiplication for adgeMatrix

Dispatches `%*%` through the amatrix backend for dense `adgeMatrix`
objects on the left-hand side, preserving GPU residency across the
operation.

## Usage

``` r
# S4 method for class 'adgeMatrix,matrix'
x %*% y

# S4 method for class 'adgeMatrix,Matrix'
x %*% y

# S4 method for class 'adgeMatrix,dgeMatrix'
x %*% y

# S4 method for class 'adgeMatrix,dgCMatrix'
x %*% y

# S4 method for class 'adgeMatrix,adgeMatrix'
x %*% y

# S4 method for class 'adgeMatrix,adgCMatrix'
x %*% y

# S4 method for class 'numeric,adgeMatrix'
x %*% y
```

## Arguments

- x:

  An `adgeMatrix`, `numeric` vector, or `matrix`.

- y:

  A matrix-like object: `matrix`, `Matrix`, `adgeMatrix`, `adgCMatrix`,
  or `ANY`.

## Value

An `adgeMatrix` (or `numeric` vector when `y` is a vector and `x` is
`adgeMatrix`), with the same backend metadata as `x`.

## Examples

``` r
A <- adgeMatrix(matrix(1:6, 2, 3))
B <- matrix(1:3, 3, 1)
A %*% B
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 2 x 1 Matrix of class "adgeMatrix"
#>      [,1]
#> [1,]   22
#> [2,]   28
```
