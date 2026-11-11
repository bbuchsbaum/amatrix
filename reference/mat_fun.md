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
#>             [,1]        [,2]        [,3]        [,4]
#> [1,]  2.62592118 -0.02566567  0.17333094 -0.05163742
#> [2,] -0.02566567  2.49365199 -0.24306464 -0.27123311
#> [3,]  0.17333094 -0.24306464  1.77703421  0.07747057
#> [4,] -0.05163742 -0.27123311  0.07747057  1.64269635
mat_log(S)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 4 x 4 Matrix of class "adgeMatrix"
#>             [,1]       [,2]       [,3]        [,4]
#> [1,]  1.92456842 -0.0150804  0.1602054 -0.05375419
#> [2,] -0.01508040  1.8004884 -0.2260995 -0.26414502
#> [3,]  0.16020537 -0.2260995  1.1264290  0.07660820
#> [4,] -0.05375419 -0.2641450  0.0766082  0.96991336
mat_pow(S, -1)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 4 x 4 Matrix of class "adgeMatrix"
#>               [,1]         [,2]        [,3]        [,4]
#> [1,]  0.1487865376 0.0003454415 -0.03638565  0.01505680
#> [2,]  0.0003454415 0.1781007915  0.05151209  0.06704050
#> [3,] -0.0363856490 0.0515120880  0.33481912 -0.02026996
#> [4,]  0.0150568021 0.0670404982 -0.02026996  0.39052370
```
