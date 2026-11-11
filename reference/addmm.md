# Scaled matrix multiply with optional bias: alpha\*(A%\*%B) + beta\*C

Scaled matrix multiply with optional bias: alpha\*(A%\*%B) + beta\*C

## Usage

``` r
addmm(A, B, C = NULL, alpha = 1, beta = 1)
```

## Arguments

- A:

  n×p `adgeMatrix` or plain matrix.

- B:

  p×k numeric matrix.

- C:

  n×k numeric matrix or `NULL` (treated as zeros).

- alpha:

  Scalar multiplier for `A%*%B` (default 1).

- beta:

  Scalar multiplier for `C` (default 1).

## Value

`adgeMatrix` if A is resident, otherwise plain matrix.

## Examples

``` r
A <- adgeMatrix(matrix(1:6, 2, 3))
B <- matrix(1:12, 3, 4)
C <- matrix(1, 2, 4)
# 2 * (A %*% B) + 0.5 * C
addmm(A, B, C, alpha = 2, beta = 0.5)
#>      [,1]  [,2]  [,3]  [,4]
#> [1,] 44.5  98.5 152.5 206.5
#> [2,] 56.5 128.5 200.5 272.5
```
