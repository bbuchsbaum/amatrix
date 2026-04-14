# Matrix functions via symmetric eigendecomposition

Apply an elementwise function to the eigenvalues of a symmetric positive
definite matrix and reconstruct the result: \\f(X) = Q \\
\mathrm{diag}(f(\lambda)) \\ Q^T\\.

## Usage

``` r
mat_sqrt(X)

mat_pow(X, p)

mat_log(X)
```

## Arguments

- X:

  Symmetric positive definite numeric matrix or `adgeMatrix` of shape
  `[p, p]`.

- p:

  Numeric scalar exponent (used by `mat_pow` only).

## Value

An `adgeMatrix` of shape `[p, p]`: the matrix function applied to `X`.

## Examples

``` r
S <- crossprod(matrix(rnorm(16), 4)) + diag(4)
mat_sqrt(S)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 4 x 4 Matrix of class "adgeMatrix"
#>             [,1]      [,2]       [,3]       [,4]
#> [1,]  2.05774913 0.9158125 -0.2264633 0.02527866
#> [2,]  0.91581253 2.7294307  0.1892996 0.25175848
#> [3,] -0.22646330 0.1892996  2.3625557 0.39730237
#> [4,]  0.02527866 0.2517585  0.3973024 1.13951865
mat_log(S)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 4 x 4 Matrix of class "adgeMatrix"
#>               [,1]      [,2]       [,3]          [,4]
#> [1,]  1.2468317429 0.8275577 -0.2517660 -0.0003514227
#> [2,]  0.8275576528 1.8375289  0.1747655  0.2657424936
#> [3,] -0.2517659518 0.1747655  1.6542121  0.4739614515
#> [4,] -0.0003514227 0.2657425  0.4739615  0.1630995714
mat_pow(S, -1)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 4 x 4 Matrix of class "adgeMatrix"
#>             [,1]        [,2]        [,3]        [,4]
#> [1,]  0.38575133 -0.21022967  0.08121069  0.01554878
#> [2,] -0.21022967  0.24570343 -0.03973318 -0.10523400
#> [3,]  0.08121069 -0.03973318  0.24069709 -0.21071536
#> [4,]  0.01554878 -0.10523400 -0.21071536  0.92834566
```
