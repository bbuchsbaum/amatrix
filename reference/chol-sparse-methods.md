# Cholesky factorization for adgCMatrix

Compute the Cholesky factorization of a sparse symmetric
positive-definite `adgCMatrix`.

## Usage

``` r
# S4 method for class 'adgCMatrix'
chol(x, ...)
```

## Arguments

- x:

  A symmetric positive-definite `adgCMatrix`.

- ...:

  Further arguments passed to
  [`Matrix::chol`](https://rdrr.io/pkg/Matrix/man/chol-methods.html).

## Value

An `adgCMatrix` or sparse Cholesky factor object.

## Examples

``` r
S <- Matrix::sparseMatrix(i = 1:3, j = 1:3, x = c(4, 9, 16))
A <- adgCMatrix(S)
chol(A)
#> An amatrix sparse matrix [cpu|policy=auto|precision=strict]
#> 3 x 3 sparse Matrix of class "adgCMatrix"
#>           
#> [1,] 2 . .
#> [2,] . 3 .
#> [3,] . . 4
```
