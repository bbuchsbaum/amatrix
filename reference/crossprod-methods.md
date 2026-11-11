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
#>            [,1]       [,2]       [,3]
#> [1,]  0.6720521 -0.3076606 -0.3906943
#> [2,] -0.3076606  2.6229662  0.5049400
#> [3,] -0.3906943  0.5049400  0.4191362
tcrossprod(A)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 4 x 4 Matrix of class "adgeMatrix"
#>            [,1]        [,2]       [,3]      [,4]
#> [1,]  0.5981518  0.12389837 -0.2905298 0.6176474
#> [2,]  0.1238984  0.09448163 -0.1219992 0.2160595
#> [3,] -0.2905298 -0.12199923  0.4280004 0.2738003
#> [4,]  0.6176474  0.21605950  0.2738003 2.5935206
```
