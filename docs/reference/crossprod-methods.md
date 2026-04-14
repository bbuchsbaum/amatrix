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
#>          [,1]     [,2]     [,3]
#> [1,] 5.391972 5.147043 2.539914
#> [2,] 5.147043 7.854458 5.169525
#> [3,] 2.539914 5.169525 4.511121
tcrossprod(A)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 4 x 4 Matrix of class "adgeMatrix"
#>          [,1]     [,2]      [,3]      [,4]
#> [1,] 3.167864 3.596732 1.8137076 3.4547569
#> [2,] 3.596732 5.630710 1.2755074 5.9257495
#> [3,] 1.813708 1.275507 1.4973132 0.7022416
#> [4,] 3.454757 5.925749 0.7022416 7.4616644
```
