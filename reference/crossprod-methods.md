# Cross-product methods for adgeMatrix

Compute \\t(x) \\\*\\ y\\ (`crossprod`) or \\x \\\*\\ t(y)\\
(`tcrossprod`) for `adgeMatrix` objects, dispatching through the amatrix
backend to preserve GPU residency.

## Usage

``` r
# S4 method for class 'adgeMatrix,ANY'
crossprod(x, y = NULL, ...)

# S4 method for class 'adgeMatrix,missing'
crossprod(x, y = NULL, ...)

# S4 method for class 'adgeMatrix,ANY'
tcrossprod(x, y = NULL, ...)

# S4 method for class 'adgeMatrix,missing'
tcrossprod(x, y = NULL, ...)
```

## Arguments

- x:

  An `adgeMatrix`.

- y:

  A matrix-like object, or `NULL` for the symmetric form \\t(x) \\\*\\
  x\\ or \\x \\\*\\ t(x)\\.

- ...:

  Further arguments passed to the underlying backend operation.

## Value

An `adgeMatrix` containing the result.

## Examples

``` r
A <- adgeMatrix(matrix(rnorm(12), 4, 3))
crossprod(A)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 3 x 3 Matrix of class "adgeMatrix"
#>           [,1]      [,2]      [,3]
#> [1,]  5.134901 -1.355748  1.840785
#> [2,] -1.355748  4.822761 -2.907902
#> [3,]  1.840785 -2.907902  3.820301
tcrossprod(A)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 4 x 4 Matrix of class "adgeMatrix"
#>          [,1]      [,2]      [,3]       [,4]
#> [1,] 6.569142  2.406681 1.6680921  1.8405693
#> [2,] 2.406681  2.905092 1.1023399 -1.0743689
#> [3,] 1.668092  1.102340 1.6400398  0.8786711
#> [4,] 1.840569 -1.074369 0.8786711  2.6636888
```
