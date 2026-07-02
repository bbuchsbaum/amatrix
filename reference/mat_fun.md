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
#>           [,1]       [,2]       [,3]       [,4]
#> [1,] 2.3878623  0.6744472  0.6235799  0.1911931
#> [2,] 0.6744472  2.3122822  1.1160848 -0.2370981
#> [3,] 0.6235799  1.1160848  1.9587087 -0.2028453
#> [4,] 0.1911931 -0.2370981 -0.2028453  1.2118061
mat_log(S)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 4 x 4 Matrix of class "adgeMatrix"
#>           [,1]       [,2]       [,3]       [,4]
#> [1,] 1.5928519  0.5036866  0.4933769  0.3026269
#> [2,] 0.5036866  1.3073505  1.0988872 -0.2712384
#> [3,] 0.4933769  1.0988872  0.9410608 -0.2388235
#> [4,] 0.3026269 -0.2712384 -0.2388235  0.3190247
mat_pow(S, -1)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 4 x 4 Matrix of class "adgeMatrix"
#>             [,1]        [,2]       [,3]       [,4]
#> [1,]  0.27007056 -0.08690596 -0.1098987 -0.1724070
#> [2,] -0.08690596  0.50234254 -0.4068124  0.1140956
#> [3,] -0.10989873 -0.40681239  0.6520366  0.1115524
#> [4,] -0.17240700  0.11409559  0.1115524  0.7886285
```
