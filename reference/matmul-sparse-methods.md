# Matrix multiplication for adgCMatrix

Dispatches `%*%` through the amatrix backend for sparse `adgCMatrix`
objects on the left-hand side, preserving GPU residency metadata across
the operation.

## Usage

``` r
# S4 method for class 'adgCMatrix,ANY'
x %*% y

# S4 method for class 'adgCMatrix,matrix'
x %*% y

# S4 method for class 'adgCMatrix,Matrix'
x %*% y

# S4 method for class 'adgCMatrix,dgeMatrix'
x %*% y

# S4 method for class 'adgCMatrix,dgCMatrix'
x %*% y

# S4 method for class 'adgCMatrix,adgeMatrix'
x %*% y

# S4 method for class 'adgCMatrix,adgCMatrix'
x %*% y
```

## Arguments

- x:

  An `adgCMatrix`.

- y:

  A matrix-like object: `matrix`, `Matrix`, `adgeMatrix`, `adgCMatrix`,
  or `ANY`.

## Value

An `adgeMatrix` or `adgCMatrix` containing the product, with backend
metadata inherited from `x`.

## Examples

``` r
sp <- as(matrix(c(1, 0, 0, 2), 2, 2), "dgCMatrix")
A  <- adgCMatrix(sp)
B  <- matrix(1:4, 2, 2)
A %*% B
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 2 x 2 Matrix of class "adgeMatrix"
#>      [,1] [,2]
#> [1,]    1    3
#> [2,]    4    8
```
