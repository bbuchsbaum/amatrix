# Cross-product methods for adgCMatrix

Compute \\t(x) \\\*\\ y\\ (`crossprod`) or \\x \\\*\\ t(y)\\
(`tcrossprod`) for sparse `adgCMatrix` objects, dispatching through the
amatrix backend.

## Usage

``` r
# S4 method for class 'adgCMatrix,missing'
crossprod(x, y = NULL, ...)

# S4 method for class 'adgCMatrix,ANY'
crossprod(x, y = NULL, ...)

# S4 method for class 'adgCMatrix,matrix'
crossprod(x, y = NULL, ...)

# S4 method for class 'adgCMatrix,Matrix'
crossprod(x, y = NULL, ...)

# S4 method for class 'adgCMatrix,dgeMatrix'
crossprod(x, y = NULL, ...)

# S4 method for class 'adgCMatrix,dgCMatrix'
crossprod(x, y = NULL, ...)

# S4 method for class 'adgCMatrix,adgeMatrix'
crossprod(x, y = NULL, ...)

# S4 method for class 'adgCMatrix,adgCMatrix'
crossprod(x, y = NULL, ...)

# S4 method for class 'adgCMatrix,missing'
tcrossprod(x, y = NULL, ...)

# S4 method for class 'adgCMatrix,ANY'
tcrossprod(x, y = NULL, ...)

# S4 method for class 'adgCMatrix,matrix'
tcrossprod(x, y = NULL, ...)

# S4 method for class 'adgCMatrix,Matrix'
tcrossprod(x, y = NULL, ...)

# S4 method for class 'adgCMatrix,dgeMatrix'
tcrossprod(x, y = NULL, ...)

# S4 method for class 'adgCMatrix,dgCMatrix'
tcrossprod(x, y = NULL, ...)

# S4 method for class 'adgCMatrix,adgeMatrix'
tcrossprod(x, y = NULL, ...)

# S4 method for class 'adgCMatrix,adgCMatrix'
tcrossprod(x, y = NULL, ...)
```

## Arguments

- x:

  An `adgCMatrix`.

- y:

  A matrix-like object or `NULL`/missing for the symmetric form.

- ...:

  Further arguments passed to the backend.

## Value

An `adgeMatrix` or `adgCMatrix` containing the result.

## Examples

``` r
sp <- as(matrix(c(1, 0, 0, 2, 1, 0), 3, 2), "dgCMatrix")
A  <- adgCMatrix(sp)
crossprod(A)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 2 x 2 Matrix of class "adgeMatrix"
#>      [,1] [,2]
#> [1,]    1    2
#> [2,]    2    5
```
