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
#>            [,1]        [,2]        [,3]       [,4]
#> [1,]  2.1941853 -0.13485827 -0.11525719  0.6386226
#> [2,] -0.1348583  1.18257693 -0.02940746 -0.4459531
#> [3,] -0.1152572 -0.02940746  1.08597129 -0.1713995
#> [4,]  0.6386226 -0.44595306 -0.17139950  2.8304459
mat_log(S)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 4 x 4 Matrix of class "adgeMatrix"
#>             [,1]        [,2]        [,3]       [,4]
#> [1,]  1.49384166 -0.09758931 -0.12197059  0.5114231
#> [2,] -0.09758931  0.25252423 -0.08830011 -0.4749923
#> [3,] -0.12197059 -0.08830011  0.14603552 -0.1830043
#> [4,]  0.51142312 -0.47499235 -0.18300425  1.9709234
mat_pow(S, -1)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 4 x 4 Matrix of class "adgeMatrix"
#>             [,1]       [,2]       [,3]        [,4]
#> [1,]  0.25222387 0.01030074 0.04577453 -0.08982148
#> [2,]  0.01030074 0.83620345 0.09764327  0.18523976
#> [3,]  0.04577453 0.09764327 0.88065489  0.07823660
#> [4,] -0.08982148 0.18523976 0.07823660  0.19649329
```
